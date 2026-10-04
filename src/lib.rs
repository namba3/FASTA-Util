#![cfg_attr(test, feature(test))]
#![feature(slice_from_ptr_range)]

pub mod amino_acid;
pub mod nucleic_acid;

pub use nucleic_acid::is_nucleic_acid_lut as is_nucleic_acid;

use core::slice;
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
    let mmap = unsafe { memmap2::Mmap::map(&file) }?;
    Ok(lines(mmap))
}

pub struct LinesInFile {
    mmap: Arc<memmap2::Mmap>,
    head: usize,
}
impl Iterator for LinesInFile {
    type Item = LineInFile;
    fn next(&mut self) -> Option<Self::Item> {
        if self.mmap.len() <= self.head {
            return None;
        }

        let line = self.mmap[self.head..]
            .split_inclusive(|byte| *byte == b'\n')
            .take(1)
            .last();

        if let Some(line) = line {
            self.head += line.len();
            let range = line.as_ptr_range();
            LineInFile {
                _mmap: Arc::clone(&self.mmap),
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
    _mmap: Arc<memmap2::Mmap>,
    slice: &'static [u8],
}
impl AsRef<[u8]> for LineInFile {
    fn as_ref<'a>(&'a self) -> &'a [u8] {
        self.slice
    }
}
fn lines(mmap: memmap2::Mmap) -> LinesInFile {
    let mmap = Arc::new(mmap);
    LinesInFile {
        mmap: mmap,
        head: 0,
    }
}

#[cfg(test)]
mod tests {
    use super::read_lines_from_file;
    use std::{
        fs::{self, File, OpenOptions},
        io::Write,
        path::PathBuf,
        sync::atomic::{AtomicUsize, Ordering},
    };

    static NEXT_TEMP_FILE_ID: AtomicUsize = AtomicUsize::new(0);

    struct TemporaryInput(PathBuf);

    impl TemporaryInput {
        fn new(contents: &[u8]) -> Self {
            loop {
                let id = NEXT_TEMP_FILE_ID.fetch_add(1, Ordering::Relaxed);
                let path = std::env::temp_dir()
                    .join(format!("fasta-util-lines-{}-{id}.tmp", std::process::id()));
                match OpenOptions::new().write(true).create_new(true).open(&path) {
                    Ok(mut file) => {
                        file.write_all(contents).unwrap();
                        return Self(path);
                    }
                    Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                    Err(error) => panic!("failed to create temporary input: {error}"),
                }
            }
        }

        fn open(&self) -> File {
            File::open(&self.0).unwrap()
        }
    }

    impl Drop for TemporaryInput {
        fn drop(&mut self) {
            let _ = fs::remove_file(&self.0);
        }
    }

    #[test]
    fn reads_each_line_including_its_newline() {
        let input = TemporaryInput::new(b">record\nACGT\nsecond line");
        let lines = read_lines_from_file(input.open())
            .unwrap()
            .map(|line| line.as_ref().to_vec())
            .collect::<Vec<_>>();

        assert_eq!(
            lines,
            [
                b">record\n".to_vec(),
                b"ACGT\n".to_vec(),
                b"second line".to_vec()
            ]
        );
    }

    #[test]
    fn does_not_yield_an_extra_line_after_a_trailing_newline() {
        let input = TemporaryInput::new(b"first\nsecond\n");
        let lines = read_lines_from_file(input.open())
            .unwrap()
            .map(|line| line.as_ref().to_vec())
            .collect::<Vec<_>>();

        assert_eq!(lines, [b"first\n".to_vec(), b"second\n".to_vec()]);
    }

    #[test]
    fn empty_file_has_no_lines() {
        let input = TemporaryInput::new(b"");

        assert_eq!(read_lines_from_file(input.open()).unwrap().count(), 0);
    }

    #[test]
    fn line_clone_keeps_mapped_bytes_alive_after_iterator_is_dropped() {
        let input = TemporaryInput::new(b"mapped line\n");
        let mut lines = read_lines_from_file(input.open()).unwrap();
        let line = lines.next().unwrap();
        let clone = line.clone();
        drop(lines);

        assert_eq!(line.as_ref(), b"mapped line\n");
        assert_eq!(clone.as_ref(), b"mapped line\n");
    }
}
