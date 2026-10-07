/// Returns whether `x` is an accepted one-letter protein sequence symbol.
///
/// The accepted set includes the 20 standard amino acids, ambiguous/rare
/// symbols (`B`, `J`, `O`, `U`, `X`, `Z`), a stop marker (`*`), and an
/// alignment gap (`-`). Lowercase symbols are accepted as well.
#[inline]
pub const fn is_amino_acid(x: u8) -> bool {
    matches!(
        x,
        b'A' | b'C'
            | b'D'
            | b'E'
            | b'F'
            | b'G'
            | b'H'
            | b'I'
            | b'K'
            | b'L'
            | b'M'
            | b'N'
            | b'P'
            | b'Q'
            | b'R'
            | b'S'
            | b'T'
            | b'V'
            | b'W'
            | b'Y'
            | b'B'
            | b'J'
            | b'O'
            | b'U'
            | b'X'
            | b'Z'
            | b'*'
            | b'-'
            | b'a'
            | b'c'
            | b'd'
            | b'e'
            | b'f'
            | b'g'
            | b'h'
            | b'i'
            | b'k'
            | b'l'
            | b'm'
            | b'n'
            | b'p'
            | b'q'
            | b'r'
            | b's'
            | b't'
            | b'v'
            | b'w'
            | b'y'
            | b'b'
            | b'j'
            | b'o'
            | b'u'
            | b'x'
            | b'z'
    )
}

#[cfg(test)]
mod tests {
    use super::is_amino_acid;

    #[test]
    fn accepts_supported_uppercase_and_lowercase_symbols() {
        for symbol in b"ACDEFGHIKLMNPQRSTVWYBJOUXZ*-acdefghiklmnpqrstvwybjouxz" {
            assert!(is_amino_acid(*symbol), "rejected {:?}", char::from(*symbol));
        }
    }

    #[test]
    fn rejects_non_protein_symbols() {
        for symbol in b"0123456789.!? " {
            assert!(
                !is_amino_acid(*symbol),
                "accepted {:?}",
                char::from(*symbol)
            );
        }
    }
}
