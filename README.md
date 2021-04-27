# numc

Here's what I did in project 4:
-Task 1, I did fairly easily
-I simply did each function without paying attention to optimization, as long as it worked, thus, none of them had any parallelization
-For each function, I wrote error handling as I went along (I would remove some later)
-PyNumberMethods and PyMethodDef were difficult for me to handle, however, and I spent a significant amount of time surfing the web to figure out how they worked
-PyNumbberMethods in particular was difficult because I didn't know for quite some time that there was a way to format it and remove the NULLs.
-deallocate was also a challenge, since I was under the impression that child matrices could have children; I ended up leaving the implementation mostly as-is, as there wasn't much to optimize there

-Task 2, too, took a significant amount of time because I didn't understand what we were supposed to do
-Eventually, I figured out that we were expected to create a module and add it as an extension
-My first iteration of setup.py was faulty though, and I would realize it until reading through piazza and learning that we were supposed to make use of the flags

-Task 3 was significantly easier than Task 2
-Task 3 was off to a bumpy start, as I had difficulty figuring how to create Matrix61c objects
-There was also a bit of a learning curve to learning how to work with PyObjects
-Like Task 1, I focused solely on implementing everything, and almost everything operated by creating an extra Matrix61c object and calling the functions created in Task 1

-Task 4 was the hardest task
-Unrolling was the first thing I implemented, as it was easy
-OpenMp was the second thing I tried to implement, and I merely added #pragma omp parallel for before any for loops
-However, this quickly ran into various issues with nested loops
-I left omp half-done and switched to SIMD
-I implemented SIMD for some functions, but I eventually chose to focus on multiplication, since I reasoned that the speed up of multiplication would affect pow and so it was the most important task to tackle
-I ended up implementing blocking first, then omp, then SIMD for multiplication
-At this point, I made very little use of #pragma omp parallel for, since I was annoyed by the issues I had encountered earlier regarding nested loops
-However, despite having the four different ways to speed up my code, I was still at only around x15 speed up for my mul and was virtually out of ideas
-I made modifications to allocate so that only one calloc had to be called, and everything about a matrix would be next to each other - including the contents of data - which should improve performance as the contents of a matrix as well as its information like rows and cols would be saved together in caches
-I came to the realization that my blocking was pretty bad - it would still exceed the cache in cases of large matrices, since it solved each element of the result matrix at a time and thus could cause conflict misses if there were too many terms to calculate 
-I tried a different way to perform mul that had to set up SIMD more often but had better blocking, in the hopes that it would be faster
-The version with worse blocking turned out to operate faster, however, so I began to revert it back
-Due to time constraints, I didn't implement blocking in the selected method in order to be able to code it faster and do more testing instead
-This also let me use #pragma omp parallel for, since the faster, short-cut version of my original code possessed fewer nested loops and didn't cause the issues that plagued me earlier
-I added SIMD and unrolling to the new multiplication function
-Funnily enough, this hastily thrown-together version did better than both of its predecessors; however, it was still slow, at only x18 speed up even with SIMD and omp implemented... which was incredibly strange, but not much I could do about it
-I had previously considered it but dismissed it as only a small increase, but since I was desperate for more ways to speed up, I implemented SIMD for the tail case... so I had a tail case and a tail case for the tail case
-My multiplication had previously also made use of fill_matrix to clean out the result matrix, as it otherwise could mess up if it was given a matrix that already had terms in it; I was able to find a way to remove this, which sped things up
-Surfing the web led me to the realization that the stack is faster than the heap, so then I proceeded to try to convert as many things as possible into variables, both to reduce the number of calculations the code had to do and to avoid memory accesses to the heap
-With the clock ticking down, I implemented the same things I did for multiplication to my other operations and switched to just testing - I was practically out of ideas for further improving mul, and was instead aiming to simply get the points for correctness