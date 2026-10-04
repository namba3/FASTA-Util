function fail(message) {
    print "cannot build FAI: " message > "/dev/stderr"
    failed = 1
    exit 1
}

function emit_record() {
    if (!have_record) {
        return
    }
    print name "\t" sequence_length "\t" sequence_offset "\t" line_bases "\t" line_width
}

{
    physical_line_width = length($0) + 1
    line = $0
    carriage_return = sub(/\r$/, "", line)
    content_length = length(line)

    if (substr(line, 1, 1) == ">") {
        emit_record()
        header = substr(line, 2)
        sub(/[[:space:]].*$/, "", header)
        if (header == "") {
            fail("empty FASTA record name")
        }
        name = header
        sequence_length = 0
        sequence_offset = byte_offset + physical_line_width
        line_bases = 0
        line_width = 0
        previous_line_bases = 0
        line_ending_width = carriage_return + 1
        have_sequence = 0
        have_record = 1
    } else {
        if (!have_record) {
            fail("sequence data appears before a FASTA header")
        }
        if (content_length == 0) {
            fail("empty sequence line")
        }
        if (carriage_return + 1 != line_ending_width) {
            fail("mixed line endings within a FASTA record")
        }
        if (!have_sequence) {
            line_bases = content_length
            line_width = content_length + line_ending_width
            have_sequence = 1
        } else if (previous_line_bases < line_bases || content_length > line_bases) {
            fail("sequence lines are not consistently wrapped")
        }
        sequence_length += content_length
        previous_line_bases = content_length
    }

    byte_offset += physical_line_width
}

END {
    if (!failed) {
        emit_record()
    }
}
