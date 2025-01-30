# Welcome to the Fundamentals of Astrodynamics Code Repository!

Folks,

It’s a great pleasure to release the codes from my book on GitHub. Over the years, I’ve received many requests for this, and now seems like the right time to do it! A few important notes:

## Code Availability and Testing

We don’t have all the routines published yet, as we’ve been testing them in each language first. This is a big task, but very fruitful as different languages have peculiarities that are beneficial in finding small errors. The [spreadsheet](Astro%20Software.xlsx) lists:
   - Routines that are essentially complete
   - Routines in progress
   - Routines that still need some work

   A few notable ones:
   - **Maneuver routines** – These are simple magnitude-only routines that are often better served by Lambert’s solution.
   - **Orbit determination (simple Batch Least Squares)** – Originally developed for SGP4DC - this needs some work to function correctly in C#. We have a beta version, but it’s not yet producing the correct outputs.
   - **Accelerations in coordinate transformations** – These are partially implemented. The C# version includes them in the overall ECI-to-ECEF transformations, and we plan to include them in all languages shortly.
   - **STK scenarios** – Many are from older versions but should still run in the latest versions of STK. We will eventually add Word documents to describe what is actually going on!

## Testing Framework

We are developing test programs but want to modify them to use the `assert` functionality available in each language, and this is still in progress. Our intent is to finalize these tests for the first set of functions before making them publicly available.

## Known Differences Between Implementations

Some routines have known discrepancies between languages, which we are working on resolving. The C# version is generally the most authoritative, followed closely by MATLAB and Python.

   Some key differences:
   - **Covariance conversions** – Not all conversions work forward and backward as indicated in our paper ([AAS 15-537](documentation/Covariance%20Transformations/AAS%2015-537%20Updated%20Analytical%20Partials%20for%20Covariance%20Transformations%20and%20OptimizationFinal.pdf)). Additionally, the retrograde factor for equinoctial elements has been changed from applying to all retrograde orbits to only those near 180-degree inclinations. This was done to minimize differences near polar orbits, but of course, it changes a lot of the results ([AAS 15-537 Rev 3](documentation/Covariance%20Transformations/AAS%2015-537%20Updated%20Analytical%20Partials%20for%20Covariance%20Transformations%20and%20OptimizationFinal%20Rev%203.pdf))
   - **Angles-only routines** – These are complex, and while some results match between languages, others do not. We are still working on these!
   - **Other discrepancies** – If you notice any issues, please let us know and open a GitHub issue! We are working on integrating these checks into our `assert`-based test framework.

## Future Plans and C++ Updates

We have included C++, but a full update is planned in the coming months. We have heard concerns about C++ memory leaks making it vulnerable to exploitation, so its long-term role in this project may have some uncertainties.

---

As always, let us know what you think, what you’d like to add, and any feedback you have. Looking forward to this new adventure with you all! 🚀

Dav