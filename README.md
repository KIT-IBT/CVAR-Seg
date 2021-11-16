# CVAR-Seg
Segmentation of S1S2 stimulation protocol measurements.

CVAR-Seg is an automated pipeline used to evaluate the S1S2 protocol. It uses the complete signal trace of each lead into account and deduces the S1 time, the number of pacing beats as well as number of pacing trains.
Finally, the S1 and S2 beats are separated and the local activation times (LATs) are evaluated.
Using the LATs the amplitude of the signal is evaluated. Additionally using the electrode locations enables space over time evaluation of CV.

An example using synthetic data is provided. Use "run_cvar_bench.m" to run the pipeline.

The complete pipeline is described in the publication "CVAR-Seg: An Automated Signal Segmentation Pipeline for Conduction Velocity and Amplitude Restitution", https://doi.org/10.3389/fphys.2021.673047
