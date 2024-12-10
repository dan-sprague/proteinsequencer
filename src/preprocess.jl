"""
    parse_proteome(path)

    Takes a `path` to a file of protein sequences and returns the set of sequences and the 99th %ile length for the set.
"""
function parse_proteome(path;minl = 500)
    seqs = FASTAReader(open(path)) do reader
		[
            (
                id = identifier(record),
                sequence = sequence(record),
                length = length(sequence(record))
            ) 
         for record in reader if length(sequence(record)) > minl
        ]
	end
        
    seqs

end 

