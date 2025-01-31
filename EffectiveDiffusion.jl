using Printf

# Constants
const R = 8.314 # Universal gas constant [J/mol·K]
const convf = 1.0E3 # Conversion factor for units
const alpha = 0.75 # Shape factor for tortuosity
const sphericity = 1.0 # Sphericity factor

# Molar masses (kg/mol or g/mol)
const Mw = [2.0, 28.0, 16.0, 44.0, 34.0, 18.0, 28.0, 32.0, 40.0, 2.0]  # H2, CO, CH4, CO2, PH3, H2O, N2, O2, Ar, He

# Species Indices (for easier reference)
const H2_index = 1
const CO_index = 2
const CH4_index = 3
const CO2_index = 4
const PH3_index = 5
const H2O_index = 6
const N2_index = 7
const O2_index = 8
const Ar_index = 9
const He_index = 10

# Molecular diameters (in Angstroms)
const sigma = [2.827, 3.690, 3.758, 3.941, 2.641, 3.798, 3.981, 3.467, 3.35, 2.57]

# Lennard-Jones energy parameters (in K)
const epsik = [59.7, 91.7, 148.6, 195.2, 809.1, 71.4, 251.5, 106.7, 143.2, 10.8]

# Calculate diffusivities and other related quantities
function general_deff(YY, Temp, radius, poro, P)
    n_species = length(Mw)
    small = 1.0E-30
    great = 1.0E+30

    # Calculating tortuosity (tao)
    tao = (1 - alpha * (1 - poro)) ^ -1

    # Calculating permeability (k)
    k = sphericity * (2.0 * radius)^2 / 150.0 * poro^3 / (1 - poro)^2

    # Knudsen diffusivity calculation
    Dk = zeros(n_species)
    for i in 1:n_species
        Dk[i] = (2.0 / 3.0) * sqrt((8.0 * R * Temp * convf) / (π * Mw[i])) * radius
    end

    # Binary diffusivity (Chapman-Enskog model)
    D_bin = zeros(n_species, n_species)
    sigma_bin = zeros(n_species, n_species)
    epsik_bin = zeros(n_species, n_species)
    Tn = zeros(n_species, n_species)
    Omega_D = zeros(n_species, n_species)
    for i in 1:n_species
        for j in 1:n_species
            sigma_bin[i, j] = (sigma[i] + sigma[j]) / 2.0
            epsik_bin[i, j] = sqrt(epsik[i] * epsik[j])
            Tn[i, j] = Temp / (epsik_bin[i, j] + small)
            Omega_D[i, j] = (1.06036 / (Tn[i, j] + small) ^ 0.15610) + (0.19300 / exp(0.47635 * Tn[i, j] + small)) + (1.03587 / exp(1.52996 * Tn[i, j] + small)) + (1.76474 / exp(3.89411 * Tn[i, j] + small))
            D_bin[i, j] = (0.001858 / 10000.0) * sqrt((Temp^3.0 * (Mw[i] + Mw[j])) / (Mw[i] * Mw[j] + small)) / (P * (sigma_bin[i, j]^2.0) * Omega_D[i, j] + small)
        end
    end

    # Molecular diffusion coefficient (Dm)
    Dm = zeros(n_species)
    for i in 1:n_species
        suma = 0.0
        for j in 1:n_species
            if i != j
                suma += YY[j] / (D_bin[i, j] + small)
            end
        end
        Dm[i] = (1.0 - YY[i]) / (suma + small)
    end

    # Effective diffusivity calculation
    Mm = sum(YY .* Mw)
    Deff = zeros(n_species)
    for i in 1:n_species
        alpha_i = 1.0 - sqrt(Mw[i] / (Mm + small))
        Deff[i] = (poro / tao) * (1.0 / ((1.0 - alpha_i * YY[i]) / (Dm[i] + small) + (1.0 / (Dk[i] + small)) + small))
    end

    return D_bin, Dm, Dk, tao, k, Deff
end

# Example usage
YY = [0.97, 0.0, 0.0, 0.0, 0.0, 0.03, 0.0, 0.0, 0.0, 0.0] # Molar fractions
Temp = 1073.0 # Temperature in K
radius = 1.07e-6 # Radius in meters
poro = 0.48 # Porosity
P = 1.0 # Pressure in atm

D_bin, Dm, Dk, tao, k, Deff = general_deff(YY, Temp, radius, poro, P)

# Display results
@show D_bin
@show Dm
@show Dk
@show tao
@show k
@show Deff
