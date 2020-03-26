function landscape_driver

%iterate over Taylor & Laurent approximation orders up to 7 in each
%direction, and output the resulting data onto the hard drive.
for i = 1:7
    for j = 1:7
        landscape_computer(i,j);
    end
end

return