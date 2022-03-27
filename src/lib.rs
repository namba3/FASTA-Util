#![feature(test)]
#![feature(once_cell)]
#![feature(slice_from_ptr_range)]

pub mod amino_acid;
pub mod nucleic_acid;

pub use nucleic_acid::is_nucleic_acid_lut as is_nucleic_acid;

use core::{ops::Deref, slice};
use std::{
    fs::File,
    io::{BufRead, BufReader},
    sync::Arc,
};

// pub enum FastaLine<'a> {
//     Header { name: &'a [u8], comment: &'a [u8] },
//     Sequence { data: &'a [u8] },
// }

pub fn read_lines_from_stdin() -> std::io::Lines<impl BufRead> {
    BufReader::new(std::io::stdin().lock()).lines()
}

pub fn read_lines_from_file(file: File) -> Result<LinesInFile, std::io::Error> {
    Ok(SharedMmap::new(file)?.lines())
}

pub struct LinesInFile {
    shared_mmap: Arc<SharedMmap>,
    head: usize,
}
impl Iterator for LinesInFile {
    type Item = LineInFile;
    fn next(&mut self) -> Option<Self::Item> {
        if self.shared_mmap.len() <= self.head {
            return None;
        }

        let line = self.shared_mmap[self.head..]
            .split_inclusive(|byte| *byte == b'\n')
            .take(1)
            .last();

        if let Some(line) = line {
            self.head += line.len();
            let range = line.as_ptr_range();
            LineInFile {
                _mmap: Arc::clone(&self.shared_mmap),
                slice: unsafe { slice::from_ptr_range::<'static, _>(range) },
            }
            .into()
        } else {
            None
        }
    }
}
#[derive(Clone)]
pub struct LineInFile {
    _mmap: Arc<SharedMmap>,
    slice: &'static [u8],
}
impl AsRef<[u8]> for LineInFile {
    fn as_ref<'a>(&'a self) -> &'a [u8] {
        self.slice
    }
}

struct SharedMmap {
    mmap: memmap::Mmap,
    _file: Arc<std::fs::File>,
}
impl Deref for SharedMmap {
    type Target = [u8];
    fn deref(&self) -> &[u8] {
        &self.mmap
    }
}
impl<'a> SharedMmap {
    fn new(file: File) -> Result<Self, std::io::Error> {
        let _file = Arc::new(file);
        let mmap = unsafe { memmap::Mmap::map(&*_file) }?;
        Ok(Self { mmap, _file })
    }
    fn lines(self) -> LinesInFile {
        let mmap = Arc::new(self);
        LinesInFile {
            shared_mmap: mmap,
            head: 0,
        }
    }
}
