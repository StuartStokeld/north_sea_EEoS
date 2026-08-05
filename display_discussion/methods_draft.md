Methods
Study system and data

We used haul-level data from the ICES North Sea International Bottom Trawl Survey (Maureaud et al., 2024), combined with reconstructed trawling effort (otter and beam, summed) from Couce et al. (2020). For each haul we calculated species richness (S), total abundance (N), and total metabolic energy (E = Σ n·mass^0.75 across length bins, normalised by the haul's minimum observed individual mass, m_min^0.75), the state variables required by the Ecological Equation of State (EEoS; Harte et al., 2022), alongside observed community biomass (B_obs).

H1: is the EEoS baseline achievable in this system?

H1 asks whether the EEoS predicts haul-level biomass from state variables alone, with no free parameters, better than chance. Our primary test compared unfitted EEoS predictions, B_pred(S,N,E), against B_obs. Because a strong test of a parameter-free theory requires bracketing its performance against both an uninformative floor and an informed ceiling, we report this alongside a benchmarking envelope: a null model permuting biomass to destroy any S/N/E–biomass relationship (floor); an unfitted productivity relationship, log(E × m_min) against log B_obs, using each haul's minimum observed individual mass as the calibration term (competitor); and a fitted log-linear model of biomass on E (ceiling). This structure evaluates whether the EEoS's advantage over a simpler productivity relationship justifies the additional theoretical machinery, not simply whether it beats chance.

As a sensitivity check on the interpretation of any absolute discrepancy between predicted and observed biomass, we tested whether the relative structure of the EEoS is recovered independent of absolute scale, following the E/B^(3/4) relationship of Harte et al. (2022, Fig. 2).
