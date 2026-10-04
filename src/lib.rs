#![cfg_attr(test, feature(test))]

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
            let slice = unsafe {
                // SAFETY: `line` is a subslice of the read-only mapping. Its pointer and length
                // describe initialized bytes in that mapping, and the cloned Arc below keeps the
                // mapping alive for as long as this slice can be accessed through `LineInFile`.
                slice::from_raw_parts::<'static, _>(line.as_ptr(), line.len())
            };
            LineInFile {
                _mmap: Arc::clone(&self.mmap),
                slice,
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

    fn collect_lines(contents: &[u8]) -> Vec<Vec<u8>> {
        let input = TemporaryInput::new(contents);
        read_lines_from_file(input.open())
            .unwrap()
            .map(|line| line.as_ref().to_vec())
            .collect()
    }

    impl Drop for TemporaryInput {
        fn drop(&mut self) {
            let _ = fs::remove_file(&self.0);
        }
    }

    #[test]
    fn reads_each_line_including_its_newline() {
        let lines = collect_lines(b">record\nACGT\nsecond line");
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
        let lines = collect_lines(b"first\nsecond\n");
        assert_eq!(lines, [b"first\n".to_vec(), b"second\n".to_vec()]);
    }

    #[test]
    fn empty_file_has_no_lines() {
        assert!(collect_lines(b"").is_empty());
    }

    #[test]
    fn preserves_empty_lines_crlf_and_a_final_line_without_newline() {
        let lines = collect_lines(b"\nfirst\r\n\nlast");

        assert_eq!(
            lines,
            [
                b"\n".to_vec(),
                b"first\r\n".to_vec(),
                b"\n".to_vec(),
                b"last".to_vec()
            ]
        );
    }

    #[test]
    fn preserves_every_byte_when_reassembling_lines() {
        let contents = (0..=u8::MAX).collect::<Vec<_>>();
        let lines = collect_lines(&contents);
        let reassembled = lines.into_iter().flatten().collect::<Vec<_>>();

        assert_eq!(reassembled, contents);
    }

    #[test]
    fn iterator_stays_exhausted_after_the_last_line() {
        let input = TemporaryInput::new(b"one line");
        let mut lines = read_lines_from_file(input.open()).unwrap();

        assert_eq!(lines.next().unwrap().as_ref(), b"one line");
        assert!(lines.next().is_none());
        assert!(lines.next().is_none());
    }

    #[test]
    fn line_clone_keeps_mapped_bytes_alive_after_iterator_is_dropped() {
        let input = TemporaryInput::new(b"mapped line\n");
        let mut lines = read_lines_from_file(input.open()).unwrap();
        let line = lines.next().unwrap();
        let clone = line.clone();
        drop(line);
        drop(lines);

        assert_eq!(clone.as_ref(), b"mapped line\n");
    }
}
