function [target_column] = increment_column(cycle_number)
alphabet = [ 'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'];
if cycle_number <= 26
    target_column = strcat(alphabet(cycle_number));
else
    m = fix(cycle_number/26);
    n = (cycle_number - m*26);
    if n == 0
        m = m - 1;
        n = 26;
    end
    target_column  = strcat(alphabet(m), alphabet(n));
end