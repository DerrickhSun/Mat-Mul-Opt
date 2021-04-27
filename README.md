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
-My first iteration of setup.py was faulty though, and I would realize it until reading through piazza and learning that we were supposed to add flags
-Task 3