pub mod amino_acid;
pub mod nucleic_acid;

pub use nucleic_acid::is_nucleic_acid_lut as is_nucleic_acid;

use std::{
    fs::File,
    io::{BufRead, BufReader},
    ops::Range,
    sync::Arc,
};

// pub enum FastaLine<'a> {
//     Header { name: &'a [u8], comment: &'a [u8] },
//     Sequence { data: &'a [u8] },
// }

pub fn read_lines_from_stdin() -> impl Iterator<Item = std::io::Result<Vec<u8>>> {
    let mut reader = BufReader::new(std::io::stdin().lock());
    let mut done = false;

    std::iter::from_fn(move || {
        if done {
            return None;
        }

        let mut line = Vec::new();
        match reader.read_until(b'\n', &mut line) {
            Ok(0) => {
                done = true;
                None
            }
            Ok(_) => {
                if line.ends_with(b"\n") {
                    line.pop();
                    if line.ends_with(b"\r") {
                        line.pop();
                    }
                }
                Some(Ok(line))
            }
            Err(error) => {
                done = true;
                Some(Err(error))
            }
        }
    })
}

/// Reads file lines through a memory map without copying the file contents.
///
/// The file must not be modified or truncated while the returned iterator or any
/// `LineInFile` yielded by it is alive. Each yielded line keeps the mapping alive,
/// even after the iterator is dropped.
pub fn read_lines_from_file(file: File) -> Result<LinesInFile, std::io::Error> {
    if file.metadata()?.len() == 0 {
        return Ok(LinesInFile {
            mmap: None,
            head: 0,
        });
    }

    // SAFETY: The returned iterator and each yielded line retain the mapping. The
    // caller must keep the backing file unchanged for those values' lifetimes,
    // as documented above.
    let mmap = unsafe { memmap2::Mmap::map(&file) }?;
    Ok(lines(mmap))
}

pub struct LinesInFile {
    mmap: Option<Arc<memmap2::Mmap>>,
    head: usize,
}
impl Iterator for LinesInFile {
    type Item = LineInFile;
    fn next(&mut self) -> Option<Self::Item> {
        let mmap = self.mmap.as_ref()?;
        let line = line_at(mmap, self.head)?;
        let start = self.head;
        self.head += line.len();

        Some(LineInFile {
            mmap: Arc::clone(mmap),
            range: start..self.head,
        })
    }
}
impl LinesInFile {
    /// Visits each line by borrowing the mapped bytes, without cloning the mapping.
    ///
    /// `line_number` is one-based. A line includes its trailing newline when present.
    pub fn try_for_each_line<E>(
        &self,
        mut visitor: impl FnMut(usize, &[u8]) -> Result<(), E>,
    ) -> Result<(), E> {
        let Some(mmap) = &self.mmap else {
            return Ok(());
        };

        let mut head = 0;
        let mut line_number = 0;
        while let Some(line) = line_at(mmap, head) {
            line_number += 1;
            visitor(line_number, line)?;
            head += line.len();
        }
        Ok(())
    }
}
#[derive(Clone)]
pub struct LineInFile {
    mmap: Arc<memmap2::Mmap>,
    range: Range<usize>,
}
impl AsRef<[u8]> for LineInFile {
    fn as_ref(&self) -> &[u8] {
        &self.mmap[self.range.clone()]
    }
}
fn lines(mmap: memmap2::Mmap) -> LinesInFile {
    let mmap = Arc::new(mmap);
    LinesInFile {
        mmap: Some(mmap),
        head: 0,
    }
}

fn line_at(mmap: &[u8], head: usize) -> Option<&[u8]> {
    if head >= mmap.len() {
        return None;
    }
    mmap[head..].split_inclusive(|byte| *byte == b'\n').next()
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
    fn borrowed_line_visitor_preserves_line_bytes_and_numbers() {
        let input = TemporaryInput::new(b"first\r\n\nlast");
        let lines = read_lines_from_file(input.open()).unwrap();
        let mut visited = Vec::new();

        lines
            .try_for_each_line(|line_number, line| {
                visited.push((line_number, line.to_vec()));
                Ok::<_, std::io::Error>(())
            })
            .unwrap();

        assert_eq!(
            visited,
            [
                (1, b"first\r\n".to_vec()),
                (2, b"\n".to_vec()),
                (3, b"last".to_vec())
            ]
        );
    }

    #[test]
    fn borrowed_line_visitor_stops_on_the_first_error() {
        let input = TemporaryInput::new(b"first\nsecond\nthird");
        let lines = read_lines_from_file(input.open()).unwrap();
        let mut visited = Vec::new();

        let error = lines
            .try_for_each_line(|line_number, line| {
                visited.push(line.to_vec());
                if line_number == 2 {
                    Err(std::io::Error::other("visitor stopped"))
                } else {
                    Ok(())
                }
            })
            .unwrap_err();

        assert_eq!(error.to_string(), "visitor stopped");
        assert_eq!(visited, [b"first\n".to_vec(), b"second\n".to_vec()]);
    }

    #[test]
    fn borrowed_line_visitor_does_not_call_back_for_an_empty_file() {
        let input = TemporaryInput::new(b"");
        let lines = read_lines_from_file(input.open()).unwrap();
        let mut visited = false;

        lines
            .try_for_each_line(|_, _| {
                visited = true;
                Ok::<_, std::io::Error>(())
            })
            .unwrap();

        assert!(!visited);
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

    #[test]
    fn lines_keep_distinct_ranges_after_iterator_is_dropped() {
        let input = TemporaryInput::new(b"same\nsame\nlast");
        let mut lines = read_lines_from_file(input.open()).unwrap();
        let first = lines.next().unwrap();
        let second = lines.next().unwrap();
        let third = lines.next().unwrap();
        drop(lines);

        assert_eq!(first.as_ref(), b"same\n");
        assert_eq!(second.as_ref(), b"same\n");
        assert_eq!(third.as_ref(), b"last");
    }
}
