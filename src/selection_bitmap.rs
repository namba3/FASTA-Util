use fasta_util::LinesInFile;
use std::io::{self, Write};

/// Compact per-record selection flags used by commands that make two passes
/// over a FASTA file.
#[derive(Default)]
pub(crate) struct SelectionBitmap {
    words: Vec<u64>,
    len: usize,
}

impl SelectionBitmap {
    pub(crate) fn push(&mut self, selected: bool) {
        let word_index = self.len / u64::BITS as usize;
        let bit_index = self.len % u64::BITS as usize;
        if word_index == self.words.len() {
            self.words.push(0);
        }
        if selected {
            self.words[word_index] |= 1_u64 << bit_index;
        }
        self.len += 1;
    }

    pub(crate) fn get(&self, index: usize) -> Option<bool> {
        if index >= self.len {
            return None;
        }
        let word_index = index / u64::BITS as usize;
        let bit_index = index % u64::BITS as usize;
        Some(self.words[word_index] & (1_u64 << bit_index) != 0)
    }
}

pub(crate) fn write_selected_records(
    lines: &LinesInFile,
    selected: &SelectionBitmap,
    writer: &mut impl Write,
) -> io::Result<()> {
    let mut record_index = None;
    lines.try_for_each_line(|_, raw_line| {
        if raw_line.first() == Some(&b'>') {
            record_index = Some(record_index.map_or(0, |index: usize| index + 1));
        }
        if record_index.is_some_and(|index| selected.get(index) == Some(true)) {
            writer.write_all(raw_line)?;
        }
        Ok(())
    })
}

#[cfg(test)]
mod tests {
    use super::SelectionBitmap;

    #[test]
    fn empty_bitmap_has_no_entries() {
        let bitmap = SelectionBitmap::default();

        assert_eq!(bitmap.get(0), None);
    }

    #[test]
    fn stores_and_reads_flags_across_word_boundaries() {
        let flags = (0..(u64::BITS as usize + 3))
            .map(|index| matches!(index, 0 | 31 | 63 | 64 | 66))
            .collect::<Vec<_>>();
        let mut bitmap = SelectionBitmap::default();
        for flag in flags.iter().copied() {
            bitmap.push(flag);
        }

        for (index, expected) in flags.iter().copied().enumerate() {
            assert_eq!(bitmap.get(index), Some(expected), "flag at index {index}");
        }
        assert_eq!(bitmap.get(flags.len()), None);
    }

    #[test]
    fn all_false_and_all_true_words_are_preserved() {
        let mut all_false = SelectionBitmap::default();
        let mut all_true = SelectionBitmap::default();
        for _ in 0..(u64::BITS as usize + 1) {
            all_false.push(false);
            all_true.push(true);
        }

        assert!((0..=u64::BITS as usize).all(|index| all_false.get(index) == Some(false)));
        assert!((0..=u64::BITS as usize).all(|index| all_true.get(index) == Some(true)));
    }
}
