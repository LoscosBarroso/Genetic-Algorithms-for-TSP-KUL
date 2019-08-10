function result = SwapSubpath(ncities, path, start, length)
i = start;
j = mod(start-1 + length - 1, ncities)+1;

while ( length > 1 ) 
   temp = path(i);
   path(i) = path(j);
   path(j) = temp;
   length = length - 2;
   if ( ++i >= ncities ) 
       i = 0;
   end
   if ( --j <  0 ) 
       j = ncities - 1;
   end
end

result = path;