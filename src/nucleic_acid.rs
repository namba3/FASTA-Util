pub const NUCLEIC_ACID_SET: &[u8; 33] = b"ACGTNUKSYMWRBDHV-acgtnuksymwrbdhv";
pub const UPPERCASE_NUCLEIC_ACID_SET: &[u8; 16] = b"ACGTNUKSYMWRBDHV";

#[inline]
pub const fn is_nucleic_acid_match(x: u8) -> bool {
    matches!(
        x,
        b'A' | b'C'
            | b'G'
            | b'T'
            | b'N'
            | b'U'
            | b'K'
            | b'S'
            | b'Y'
            | b'M'
            | b'W'
            | b'R'
            | b'B'
            | b'D'
            | b'H'
            | b'V'
            | b'-'
            | b'a'
            | b'c'
            | b'g'
            | b't'
            | b'n'
            | b'u'
            | b'k'
            | b's'
            | b'y'
            | b'm'
            | b'w'
            | b'r'
            | b'b'
            | b'd'
            | b'h'
            | b'v'
    )
}

#[inline]
pub fn is_nucleic_acid_iter(x: u8) -> bool {
    NUCLEIC_ACID_SET.contains(&x)
}

#[inline]
pub const fn is_nucleic_acid_lut(x: u8) -> bool {
    const LUT: [bool; 256] = {
        let mut v = [false; 256];
        v[b'A' as usize] = true;
        v[b'C' as usize] = true;
        v[b'G' as usize] = true;
        v[b'T' as usize] = true;
        v[b'N' as usize] = true;
        v[b'U' as usize] = true;

        v[b'K' as usize] = true;
        v[b'S' as usize] = true;
        v[b'Y' as usize] = true;
        v[b'M' as usize] = true;
        v[b'W' as usize] = true;
        v[b'R' as usize] = true;
        v[b'B' as usize] = true;
        v[b'D' as usize] = true;
        v[b'H' as usize] = true;
        v[b'V' as usize] = true;
        v[b'-' as usize] = true;
        v[b'a' as usize] = true;
        v[b'c' as usize] = true;
        v[b'g' as usize] = true;
        v[b't' as usize] = true;
        v[b'n' as usize] = true;
        v[b'u' as usize] = true;
        v[b'k' as usize] = true;
        v[b's' as usize] = true;
        v[b'y' as usize] = true;
        v[b'm' as usize] = true;
        v[b'w' as usize] = true;
        v[b'r' as usize] = true;
        v[b'b' as usize] = true;
        v[b'd' as usize] = true;
        v[b'h' as usize] = true;
        v[b'v' as usize] = true;

        v
    };

    LUT[x as usize]
}

#[cfg(test)]
mod tests {
    macro_rules! test_is_nucleobase {
        ($name:ident) => {
            mod $name {
                use super::super::NUCLEIC_ACID_SET;
                use super::super::$name;

                #[test]
                fn accept_valid_elems() {
                    for x in NUCLEIC_ACID_SET {
                        assert!($name(*x));
                    }
                }

                #[test]
                fn reject_invalid_elems() {
                    for x in (0..=255u8).filter(|x| !NUCLEIC_ACID_SET.contains(x)) {
                        assert!(!$name(x));
                    }
                }
            }
        };
    }

    test_is_nucleobase!(is_nucleic_acid_match);
    test_is_nucleobase!(is_nucleic_acid_iter);
    test_is_nucleobase!(is_nucleic_acid_lut);
}
