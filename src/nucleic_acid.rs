#![allow(dead_code)]

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[repr(u8)]
enum NucleicAcid {
    /// Adenosine
    A = b'A',
    /// Cytidine
    C = b'C',
    /// Guanine
    G = b'G',
    /// Thymidine
    T = b'T',
    /// A/G/C/T
    N = b'N',
    /// Uridine
    U = b'U',
    /// G/T (keto)
    K = b'K',
    /// G/C (strong)
    S = b'S',
    /// T/C (pyrimidine)
    Y = b'Y',
    /// A/C (amino)
    M = b'M',
    /// A/T (weak)
    W = b'W',
    /// G/A (purine)
    R = b'R',
    /// G/T/C
    B = b'B',
    /// G/A/T
    D = b'D',
    /// A/C/T
    H = b'H',
    /// G/C/A
    V = b'V',
    /// Gap of indeterminate length
    GAP = b'-',
}
impl PartialEq<u8> for NucleicAcid {
    fn eq(&self, rhs: &u8) -> bool {
        *self as u8 == *rhs
    }
}
impl Into<u8> for NucleicAcid {
    fn into(self) -> u8 {
        self as u8
    }
}

pub const NUCLEIC_ACID_SET: &[u8; 17] = b"ACGTNUKSYMWRBDHV-";

#[inline]
pub const fn is_nucleic_acid_match(x: u8) -> bool {
    match x {
        b'A' | b'C' | b'G' | b'T' | b'N' | b'U' => true,
        b'K' | b'S' | b'Y' | b'M' | b'W' | b'R' | b'B' | b'D' | b'H' | b'V' | b'-' => true,
        _ => false,
    }
}

#[inline]
pub fn is_nucleic_acid_iter(x: u8) -> bool {
    NUCLEIC_ACID_SET.iter().find(|e| **e == x).is_some()
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

        v
    };

    LUT[x as usize]
}

#[cfg(test)]
mod tests {
    macro_rules! test_is_nucleobase {
        ($name:ident) => {
            mod $name {
                use super::super::$name;
                use super::super::NUCLEIC_ACID_SET;

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

#[cfg(test)]
mod benchs {
    extern crate test;
    use super::NUCLEIC_ACID_SET;
    use rand::Rng;
    use std::lazy::SyncLazy;
    use test::Bencher;

    static SEQ: SyncLazy<Vec<u8>> = SyncLazy::new(|| {
        const SIZE: usize = 10000;

        let mut rand = rand::thread_rng();

        let mut v = Vec::with_capacity(SIZE);
        for _ in 0..SIZE {
            let idx = rand.gen::<usize>() % NUCLEIC_ACID_SET.len();
            v.push(NUCLEIC_ACID_SET[idx]);
        }

        v
    });

    macro_rules! bench_is_nucleobase {
        ($name:ident) => {
            #[bench]
            fn $name(b: &mut Bencher) {
                let seq = &SEQ[..];
                b.iter(|| seq.iter().copied().filter(|x| super::$name(*x)).count())
            }
        };
    }

    bench_is_nucleobase!(is_nucleic_acid_match);
    bench_is_nucleobase!(is_nucleic_acid_iter);
    bench_is_nucleobase!(is_nucleic_acid_lut);
}
