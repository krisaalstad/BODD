# BODD
Bayesian Optimization for a Degree Day model (BODD)
(Coincidence: BODD also means "lived" in Norwegian)

This is a minimal working example of Bayesian Optimization (BO) for calibrating cryospheric models in pure Matlab code originally created during the Machine Learning session at the It's The Cryogrid Hackathon (ITCH) in 2024. I put it on github in case it is useful now more than a year later, and as a lesson for the future I probably should have done this sooner. Hopefully it is helpful for gaining an understanding of how BO works and maybe even someone wants to remix this so as to use it with more complex models like CryoGrid. The benefit of using Matlab is you would not need a wrapper so it could be done as a part of `inner loops'.

The method is implemented in a bare bones manner here but the overall recipe is simple: train an emulator (surrogate model) typically in the form of a Gaussian Process (GP) to mimic the parameter response of your model in terms of some objective (pick your favorite, from Mean Squared Error to something more exotic) which can then be used to pick new locations (ideally based on Bayesian Decision Theory, hence the `Bayesian') in parameter space to run your full model. This can converge in very few steps, is fully automatable, and generally much more efficient than other brute force techniques (grid search) or hand tuning. It is one of the standard methods for tuning hyperparameters in deep learning in a field known as autoML.

For more on BO (and GPs) I highly recommend the free book by Garnett:
[Bayesian Optimizatin](https://bayesoptbook.com/)
and for an introduction to Gaussian Processes ("infinite width neural networks" if you want to show off...) check out the classic GP bible by Rasmussen and Williams which is also freely available here:
[Gaussian Processes for Machine Learning](https://gaussianprocess.org/gpml/)

Don't hesitate to get in touch if you have questions, want to work on something related to Bayesian methods together, or spot mistakes in the code which there are surely many of.

Cheers,
Kris 




