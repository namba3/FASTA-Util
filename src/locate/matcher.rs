use std::{collections::VecDeque, io};

pub(super) const MAX_BIT_PARALLEL_MISMATCHES: usize = 3;

pub(super) fn parse_pattern(pattern: &[u8]) -> io::Result<(Vec<u8>, Vec<u8>)> {
    if pattern.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "motif cannot be empty",
        ));
    }
    let masks = pattern
        .iter()
        .copied()
        .map(|byte| {
            iupac_mask(byte).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "invalid IUPAC nucleotide symbol '{}' in motif",
                        char::from(byte)
                    ),
                )
            })
        })
        .collect::<io::Result<Vec<_>>>()?;
    let reverse = masks.iter().rev().copied().map(complement_mask).collect();
    Ok((masks, reverse))
}

pub(super) struct BitParallelMatcher {
    forward: ShiftAndState,
    reverse: Option<ShiftAndState>,
}

impl BitParallelMatcher {
    pub(super) fn new(forward_pattern: &[u8], reverse_pattern: &[u8], max_mismatch: usize) -> Self {
        Self {
            forward: ShiftAndState::new(forward_pattern, max_mismatch),
            reverse: (forward_pattern != reverse_pattern)
                .then(|| ShiftAndState::new(reverse_pattern, max_mismatch)),
        }
    }

    pub(super) fn reset(&mut self) {
        self.forward.reset();
        if let Some(reverse) = &mut self.reverse {
            reverse.reset();
        }
    }

    pub(super) fn advance(&mut self, sequence_mask: u8) -> (bool, bool) {
        let forward = self.forward.advance(sequence_mask);
        let reverse = self
            .reverse
            .as_mut()
            .map_or(forward, |matcher| matcher.advance(sequence_mask));
        (forward, reverse)
    }
}

enum ShiftAndState {
    Single(SingleWordShiftAndState),
    MultiWord(Box<MultiWordShiftAndState>),
}

struct SingleWordShiftAndState {
    matching_positions: [u64; 17],
    exact_state: u64,
    one_mismatch_state: u64,
    final_position: u64,
    max_mismatch: usize,
}

impl ShiftAndState {
    fn new(pattern: &[u8], max_mismatch: usize) -> Self {
        if pattern.len() > u64::BITS as usize || max_mismatch > 1 {
            return Self::MultiWord(Box::new(MultiWordShiftAndState::new(pattern, max_mismatch)));
        }

        let mut matching_positions = [0u64; 17];
        for (index, pattern_mask) in pattern.iter().copied().enumerate() {
            let position = 1u64 << index;
            for (sequence_mask, matches) in matching_positions.iter_mut().enumerate().skip(1) {
                if pattern_mask & sequence_mask as u8 != 0 {
                    *matches |= position;
                }
            }
        }
        Self::Single(SingleWordShiftAndState {
            matching_positions,
            exact_state: 0,
            one_mismatch_state: 0,
            final_position: 1u64 << (pattern.len() - 1),
            max_mismatch,
        })
    }

    fn reset(&mut self) {
        match self {
            Self::Single(state) => state.reset(),
            Self::MultiWord(state) => state.reset(),
        }
    }

    fn advance(&mut self, sequence_mask: u8) -> bool {
        match self {
            Self::Single(state) => state.advance(sequence_mask),
            Self::MultiWord(state) => state.advance(sequence_mask),
        }
    }
}

impl SingleWordShiftAndState {
    fn reset(&mut self) {
        self.exact_state = 0;
        self.one_mismatch_state = 0;
    }

    fn advance(&mut self, sequence_mask: u8) -> bool {
        let matching_positions = self.matching_positions[sequence_mask as usize];
        let previous_exact = self.exact_state;
        self.exact_state = ((previous_exact << 1) | 1) & matching_positions;
        if self.max_mismatch == 0 {
            return self.exact_state & self.final_position != 0;
        }

        let previous_one_mismatch = self.one_mismatch_state;
        self.one_mismatch_state =
            (((previous_one_mismatch << 1) | 1) & matching_positions) | ((previous_exact << 1) | 1);
        self.one_mismatch_state & self.final_position != 0
    }
}

struct MultiWordShiftAndState {
    matching_positions: [Vec<u64>; 17],
    mode: MultiWordMode,
    max_mismatch: usize,
    final_word: usize,
    final_position: u64,
}

enum MultiWordMode {
    UpToOne {
        exact_state: Vec<u64>,
        one_mismatch_state: Vec<u64>,
    },
    Multiple {
        states_by_word: Vec<u64>,
    },
}

impl MultiWordShiftAndState {
    fn new(pattern: &[u8], max_mismatch: usize) -> Self {
        let word_count = pattern.len().div_ceil(u64::BITS as usize);
        let mut matching_positions = std::array::from_fn(|_| vec![0u64; word_count]);
        for (index, pattern_mask) in pattern.iter().copied().enumerate() {
            let word = index / u64::BITS as usize;
            let position = 1u64 << (index % u64::BITS as usize);
            for (sequence_mask, words) in matching_positions.iter_mut().enumerate().skip(1) {
                if pattern_mask & sequence_mask as u8 != 0 {
                    words[word] |= position;
                }
            }
        }
        let mode = if max_mismatch <= 1 {
            MultiWordMode::UpToOne {
                exact_state: vec![0; word_count],
                one_mismatch_state: vec![0; word_count],
            }
        } else {
            MultiWordMode::Multiple {
                states_by_word: vec![0; word_count * (max_mismatch + 1)],
            }
        };
        Self {
            matching_positions,
            mode,
            max_mismatch,
            final_word: (pattern.len() - 1) / u64::BITS as usize,
            final_position: 1u64 << ((pattern.len() - 1) % u64::BITS as usize),
        }
    }

    fn reset(&mut self) {
        match &mut self.mode {
            MultiWordMode::UpToOne {
                exact_state,
                one_mismatch_state,
            } => {
                exact_state.fill(0);
                one_mismatch_state.fill(0);
            }
            MultiWordMode::Multiple { states_by_word } => states_by_word.fill(0),
        }
    }

    fn advance(&mut self, sequence_mask: u8) -> bool {
        let matching_positions = &self.matching_positions[sequence_mask as usize];
        let max_mismatch = self.max_mismatch;
        let final_word = self.final_word;
        let final_position = self.final_position;
        match &mut self.mode {
            MultiWordMode::UpToOne {
                exact_state,
                one_mismatch_state,
            } => {
                let mut exact_carry = 1;
                let mut one_mismatch_carry = 1;
                for (word, matching_positions) in matching_positions.iter().copied().enumerate() {
                    let previous_exact = exact_state[word];
                    let shifted_exact = (previous_exact << 1) | exact_carry;
                    exact_state[word] = shifted_exact & matching_positions;
                    exact_carry = previous_exact >> (u64::BITS - 1);
                    if max_mismatch == 1 {
                        let previous_one_mismatch = one_mismatch_state[word];
                        one_mismatch_state[word] = (((previous_one_mismatch << 1)
                            | one_mismatch_carry)
                            & matching_positions)
                            | shifted_exact;
                        one_mismatch_carry = previous_one_mismatch >> (u64::BITS - 1);
                    }
                }
                let state = if max_mismatch == 0 {
                    exact_state
                } else {
                    one_mismatch_state
                };
                state[final_word] & final_position != 0
            }
            MultiWordMode::Multiple { states_by_word } => {
                let mut carries = [1u64; MAX_BIT_PARALLEL_MISMATCHES + 1];
                let mut next_carries = [0u64; MAX_BIT_PARALLEL_MISMATCHES + 1];
                let state_count = max_mismatch + 1;
                for (word, matching_positions) in matching_positions.iter().copied().enumerate() {
                    let mut previous_states = [0u64; MAX_BIT_PARALLEL_MISMATCHES + 1];
                    let state_offset = word * state_count;
                    for (error_count, previous_state) in
                        previous_states.iter_mut().enumerate().take(state_count)
                    {
                        *previous_state = states_by_word[state_offset + error_count];
                    }
                    for error_count in 0..=max_mismatch {
                        let previous_state = previous_states[error_count];
                        let shifted_with_match = (previous_state << 1) | carries[error_count];
                        let next_state = if error_count == 0 {
                            shifted_with_match & matching_positions
                        } else {
                            (shifted_with_match & matching_positions)
                                | ((previous_states[error_count - 1] << 1)
                                    | carries[error_count - 1])
                        };
                        states_by_word[state_offset + error_count] = next_state;
                        next_carries[error_count] = previous_state >> (u64::BITS - 1);
                    }
                    carries = next_carries;
                }
                states_by_word[final_word * state_count + max_mismatch] & final_position != 0
            }
        }
    }
}

pub(super) fn window_matches(
    window: &VecDeque<u8>,
    forward_pattern: &[u8],
    reverse_pattern: &[u8],
    max_mismatch: usize,
) -> (bool, bool) {
    if forward_pattern == reverse_pattern {
        let mut mismatches = 0;
        for (sequence_mask, pattern_mask) in window.iter().copied().zip(forward_pattern.iter()) {
            mismatches += usize::from(sequence_mask & pattern_mask == 0);
            if mismatches > max_mismatch {
                return (false, false);
            }
        }
        let matched = mismatches <= max_mismatch;
        return (matched, matched);
    }

    let mut forward_mismatches = 0;
    let mut reverse_mismatches = 0;
    for ((sequence_mask, forward_mask), reverse_mask) in window
        .iter()
        .copied()
        .zip(forward_pattern.iter().copied())
        .zip(reverse_pattern.iter().copied())
    {
        forward_mismatches += usize::from(sequence_mask & forward_mask == 0);
        reverse_mismatches += usize::from(sequence_mask & reverse_mask == 0);
        if forward_mismatches > max_mismatch && reverse_mismatches > max_mismatch {
            break;
        }
    }
    (
        forward_mismatches <= max_mismatch,
        reverse_mismatches <= max_mismatch,
    )
}

pub(super) fn iupac_mask(byte: u8) -> Option<u8> {
    Some(match byte.to_ascii_uppercase() {
        b'A' => 0b0001,
        b'C' => 0b0010,
        b'G' => 0b0100,
        b'T' | b'U' => 0b1000,
        b'R' => 0b0101,
        b'Y' => 0b1010,
        b'S' => 0b0110,
        b'W' => 0b1001,
        b'K' => 0b1100,
        b'M' => 0b0011,
        b'B' => 0b1110,
        b'D' => 0b1101,
        b'H' => 0b1011,
        b'V' => 0b0111,
        b'N' => 0b1111,
        b'-' => 0b1_0000,
        _ => return None,
    })
}

fn complement_mask(mask: u8) -> u8 {
    ((mask & 0b0001) << 3)
        | ((mask & 0b0010) << 1)
        | ((mask & 0b0100) >> 1)
        | ((mask & 0b1000) >> 3)
        | (mask & 0b1_0000)
}

#[cfg(test)]
mod tests {
    use super::{BitParallelMatcher, complement_mask, iupac_mask, parse_pattern, window_matches};
    use crate::is_nucleic_acid;
    use std::collections::VecDeque;

    #[test]
    fn motif_parser_rejects_empty_and_invalid_patterns() {
        let empty = parse_pattern(b"").unwrap_err();
        let invalid = parse_pattern(b"ACZ").unwrap_err();

        assert_eq!(empty.kind(), std::io::ErrorKind::InvalidInput);
        assert_eq!(invalid.kind(), std::io::ErrorKind::InvalidInput);
        assert!(
            invalid
                .to_string()
                .contains("invalid IUPAC nucleotide symbol 'Z'")
        );
    }

    #[test]
    fn maps_standard_iupac_symbols_to_base_sets() {
        let symbols = b"ACGTURYSWKMBDHVN-";
        let masks = symbols
            .iter()
            .map(|symbol| iupac_mask(*symbol).unwrap())
            .collect::<Vec<_>>();

        assert_eq!(
            masks,
            [
                0b0001, 0b0010, 0b0100, 0b1000, 0b1000, 0b0101, 0b1010, 0b0110, 0b1001, 0b1100,
                0b0011, 0b1110, 0b1101, 0b1011, 0b0111, 0b1111, 0b1_0000,
            ]
        );
    }

    #[test]
    fn iupac_mask_coverage_matches_the_nucleotide_validator() {
        for byte in 0..=u8::MAX {
            assert_eq!(
                iupac_mask(byte).is_some(),
                is_nucleic_acid(byte),
                "byte {byte}"
            );
        }
    }

    #[test]
    fn complementing_an_iupac_mask_twice_returns_the_original_set() {
        for symbol in b"ACGTURYSWKMBDHVN-" {
            let mask = iupac_mask(*symbol).unwrap();
            assert_eq!(complement_mask(complement_mask(mask)), mask);
        }
    }

    #[test]
    fn palindromic_patterns_match_both_strands_with_one_mismatch_count() {
        let pattern = [0b0001, 0b0010, 0b0010, 0b0001];
        let exact_window = VecDeque::from(pattern);
        let one_mismatch_window = VecDeque::from([0b0001, 0b0010, 0b0100, 0b0001]);
        let two_mismatch_window = VecDeque::from([0b0100, 0b0010, 0b0100, 0b0001]);

        assert_eq!(
            window_matches(&exact_window, &pattern, &pattern, 0),
            (true, true)
        );
        assert_eq!(
            window_matches(&one_mismatch_window, &pattern, &pattern, 1),
            (true, true)
        );
        assert_eq!(
            window_matches(&two_mismatch_window, &pattern, &pattern, 1),
            (false, false)
        );
    }

    #[test]
    fn bit_parallel_matching_agrees_with_window_matching_for_iupac_masks() {
        let mut patterns = vec![
            vec![0b0001],
            vec![0b0001, 0b1000],
            vec![0b0101, 0b1010],
            vec![0b0001, 0b0010, 0b0100],
            vec![0b0101, 0b1111, 0b1010, 0b1_0000],
            vec![0b0001; 64],
        ];
        patterns.extend([63, 64, 65, 127, 128, 129].map(|length| {
            (0..length)
                .map(|index| [1, 2, 4, 8, 5, 10, 15, 16, 3, 12][index % 10])
                .collect()
        }));
        let sequence = (0..257)
            .map(|index| [1, 2, 4, 8, 5, 10, 15, 16, 3, 12][index % 10])
            .collect::<Vec<_>>();

        for pattern in patterns {
            let reverse_pattern = pattern
                .iter()
                .rev()
                .copied()
                .map(complement_mask)
                .collect::<Vec<_>>();
            for max_mismatch in 0..=3.min(pattern.len()) {
                let mut matcher = BitParallelMatcher::new(&pattern, &reverse_pattern, max_mismatch);
                let mut window = VecDeque::with_capacity(pattern.len());

                for (position, mask) in sequence.iter().copied().enumerate() {
                    let actual = matcher.advance(mask);
                    window.push_back(mask);
                    if window.len() > pattern.len() {
                        window.pop_front();
                    }
                    let expected = if window.len() == pattern.len() {
                        window_matches(&window, &pattern, &reverse_pattern, max_mismatch)
                    } else {
                        (false, false)
                    };
                    assert_eq!(
                        actual,
                        expected,
                        "pattern length {}, position {position}, mismatch {max_mismatch}",
                        pattern.len()
                    );
                }
            }
        }
    }

    #[test]
    fn multiword_multiple_mismatch_state_resets_between_records() {
        let pattern = [0b0001; 65];
        let reverse_pattern = [0b1000; 65];
        let mut matcher = BitParallelMatcher::new(&pattern, &reverse_pattern, 3);
        let mut result = (false, false);

        for _ in 0..pattern.len() {
            result = matcher.advance(0b0001);
        }
        assert_eq!(result, (true, false));

        matcher.reset();
        for _ in 0..pattern.len() {
            result = matcher.advance(0b0010);
        }
        assert_eq!(result, (false, false));
    }
}
