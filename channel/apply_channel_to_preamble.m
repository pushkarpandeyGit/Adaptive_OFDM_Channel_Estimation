function receivedPreamble = apply_channel_to_preamble(preamble, h)

preamble = preamble(:).';
h = h(:).';

if isempty(preamble)
    error('Preamble cannot be empty.');
end

if isempty(h)
    error('Channel impulse response cannot be empty.');
end

receivedPreamble = apply_channel(preamble, h);

end