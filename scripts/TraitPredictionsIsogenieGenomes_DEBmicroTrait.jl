using DEBmicroTrait
using CSV, DataFrames
using Roots
using Statistics
using Plots

df     = CSV.read("/Users/glmarschmann/Data/Zhen/IsogenieGenomes.ecosysguilds_relab.csv", DataFrame)

################################################################################################################################################################

# 1 & 2 Aerobic heterotroph
df_aerob_hetero                 = df[(df.heterotroph_withETC.==1) , :]
Genome_size_aerob_hetero        = convert(Array{Float64}, df_aerob_hetero.bp)
V_cell_aerob_hetero             = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_aerob_hetero)
Min_gen_time_aerob_hetero       = df_aerob_hetero.minimum_generation_time_hours
rrn_copies_aerob_hetero         = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_aerob_hetero)
Gram_stain_aerob_hetero         = repeat(["+"], size(V_cell_aerob_hetero,1))

################################################################################
# Glucose
N_C  = [6]
y_DE = [3.09]./N_C
y_EM = ones(size(V_cell_aerob_hetero,1))
ρ_ps_aerob_hetero_glucose = zeros(size(V_cell_aerob_hetero,1))
r_ms = zeros(size(V_cell_aerob_hetero,1))
for i in 1:size(V_cell_aerob_hetero,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_aerob_hetero[i]], [Min_gen_time_aerob_hetero[i]], [Gram_stain_aerob_hetero[i]], [rrn_copies_aerob_hetero[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_aerob_hetero_glucose[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_aerob_hetero[i]], [Min_gen_time_aerob_hetero[i]], [Gram_stain_aerob_hetero[i]], [rrn_copies_aerob_hetero[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV                           = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_aerob_hetero)
r_growth                       = @. (1-1/y_EV)*r_ms
r_main                         = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_aerob_hetero, Min_gen_time_aerob_hetero, Gram_stain_aerob_hetero)
N_SB_aerob_hetero_glucose      = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_aerob_hetero[i]], ρ_ps_aerob_hetero_glucose[i].*ones(1,1), [Min_gen_time_aerob_hetero[i]], [Gram_stain_aerob_hetero[i]])[1] for i in 1:size(V_cell_aerob_hetero,1)]
k2p                            = 180.0*60^2
r_assim                        = @. y_DE*N_C*k2p*N_SB_aerob_hetero_glucose
r_resp_aerob_hetero_glucose    = r_growth .+ r_main .+r_assim
median(r_resp_aerob_hetero_glucose)
std(r_resp_aerob_hetero_glucose)
r_resp_aerob_hetero_glucose_scaled = r_resp_aerob_hetero_glucose.*df_aerob_hetero.medab
median(r_resp_aerob_hetero_glucose_scaled)
std(r_resp_aerob_hetero_glucose_scaled)
histogram(r_resp_aerob_hetero_glucose)
xlabel!("Max. respiration rate [1/h]")
histogram(r_resp_aerob_hetero_glucose_scaled)
xlabel!("Max. respiration rate scaled [1/h]")
histogram(df_aerob_hetero.medab)
xlabel!("DNA abundance")
# half-saturation constant
Molecular_weight       = [180.16]
D_S                    = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_aerob_hetero       = DEBmicroTrait.specific_reference_affinity(V_cell_aerob_hetero, reshape(ρ_ps_aerob_hetero_glucose*100, 304,1), D_S)
median(K_D_aerob_hetero.*N_C.*12.011)
std(K_D_aerob_hetero.*N_C*12.011)
K_D_aerob_hetero_glucose = K_D_aerob_hetero.*N_C.*12.011
################################################################################
# Acetate
N_C  = [2]
y_DE = [0.93]./N_C
y_EM = ones(size(V_cell_aerob_hetero,1))
ρ_ps_aerob_hetero_acetate = zeros(size(V_cell_aerob_hetero,1))
r_ms = zeros(size(V_cell_aerob_hetero,1))
for i in 1:size(V_cell_aerob_hetero,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_aerob_hetero[i]], [Min_gen_time_aerob_hetero[i]], [Gram_stain_aerob_hetero[i]], [rrn_copies_aerob_hetero[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_aerob_hetero_acetate[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_aerob_hetero[i]], [Min_gen_time_aerob_hetero[i]], [Gram_stain_aerob_hetero[i]], [rrn_copies_aerob_hetero[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV                           = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_aerob_hetero)
r_growth                       = @. (1-1/y_EV)*r_ms
r_main                         = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_aerob_hetero, Min_gen_time_aerob_hetero, Gram_stain_aerob_hetero)
N_SB_aerob_hetero_acetate      = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_aerob_hetero[i]], ρ_ps_aerob_hetero_acetate[i].*ones(1,1), [Min_gen_time_aerob_hetero[i]], [Gram_stain_aerob_hetero[i]])[1] for i in 1:size(V_cell_aerob_hetero,1)]
k2p                            = 180.0*60^2
r_assim                        = @. y_DE*N_C*k2p*N_SB_aerob_hetero_acetate
r_resp_aerob_hetero_acetate    = r_growth .+ r_main .+r_assim
median(r_resp_aerob_hetero_acetate)
std(r_resp_aerob_hetero_acetate)
# half-saturation constant
Molecular_weight       = [59.04]
D_S                    = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_aerob_hetero       = DEBmicroTrait.specific_reference_affinity(V_cell_aerob_hetero, reshape(ρ_ps_aerob_hetero_acetate*100, 304,1), D_S)
median(K_D_aerob_hetero.*N_C.*12.011)
std(K_D_aerob_hetero.*N_C*12.011)
K_D_aerob_hetero_acetate = K_D_aerob_hetero.*N_C.*12.011

df_out = DataFrame()
df_out.genome = vcat(df_aerob_hetero.genome, df_aerob_hetero.genome)
df_out.domaine = vcat(df_aerob_hetero.domain, df_aerob_hetero.domain)
df_out.phylum = vcat(df_aerob_hetero.phylum, df_aerob_hetero.phylum)
df_out.habitat_binnedfrom = vcat(df_aerob_hetero.habitat_binnedfrom, df_aerob_hetero.habitat_binnedfrom)
df_out.depth_binnedfrom = vcat(df_aerob_hetero.depth_binnedfrom, df_aerob_hetero.depth_binnedfrom)
df_out.functional_group = repeat(["Aerobic heterotroph"], 2*304)
df_out.substrate = vcat(repeat(["Glucose"], 304), repeat(["Acetate"], 304))
df_out.Rmax = vcat(r_resp_aerob_hetero_glucose, r_resp_aerob_hetero_acetate)
df_out.K = vcat(K_D_aerob_hetero_glucose[:,1], K_D_aerob_hetero_acetate[:,1])
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/emerge/files/TraitPredictionsIsogenieGenomes_1and2AerobicHeterotophs.csv", df_out)
################################################################################################################################################################

################################################################################################################################################################
# 6 Methanotroph

df_mo = df[(df.methane_oxidation.==1) , :]
Genome_size_mo = convert(Array{Float64}, df_mo.bp)
V_cell_mo = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_mo)
Min_gen_time_mo = df_mo.minimum_generation_time_hours
rrn_copies_mo = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_mo)
Gram_stain_mo = repeat(["+"], size(V_cell_mo,1))

N_C  = [1]
y_DE = [0.88]./N_C
y_EM = ones(size(V_cell_mo,1))
ρ_ps_mo = zeros(size(V_cell_mo,1))
r_ms = zeros(size(V_cell_mo,1))
for i in 1:size(V_cell_mo,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_mo[i]], [Min_gen_time_mo[i]], [Gram_stain_mo[i]], [rrn_copies_mo[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_mo[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_mo[i]], [Min_gen_time_mo[i]], [Gram_stain_mo[i]], [rrn_copies_mo[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV                           = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_mo)
r_growth                       = @. (1-1/y_EV)*r_ms
r_main                         = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_mo, Min_gen_time_mo, Gram_stain_mo)
N_SB_mo                        = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_mo[i]], ρ_ps_mo[i].*ones(1,1), [Min_gen_time_mo[i]], [Gram_stain_mo[i]])[1] for i in 1:size(V_cell_mo,1)]
k2p                            = 180.0*60^2
r_assim                        = @. y_DE*N_C*k2p*N_SB_mo
r_resp_mo     = r_growth .+ r_main .+r_assim
median(r_resp_mo)
std(r_resp_mo)
# half-saturation constant
Molecular_weight       = [16.043]
D_S                    = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_mo       = DEBmicroTrait.specific_reference_affinity(V_cell_mo, reshape(ρ_ps_mo/100, 9,1), D_S)
median(K_D_mo.*N_C.*12.011)
std(K_D_mo.*N_C*12.011)
K_D_mo= K_D_mo.*N_C.*12.011

df_out = DataFrame()
df_out.genome = vcat(df_mo.genome)
df_out.domaine = vcat(df_mo.domain)
df_out.phylum = vcat(df_mo.phylum)
df_out.habitat_binnedfrom = vcat(df_mo.habitat_binnedfrom)
df_out.depth_binnedfrom = vcat(df_mo.depth_binnedfrom)
df_out.functional_group = repeat(["Methanotroph"], 9)
df_out.substrate = vcat(repeat(["Methane"], 9))
df_out.Rmax = vcat(r_resp_mo)
df_out.K = vcat(K_D_mo[:,1])
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/emerge/files/TraitPredictionsIsogenieGenomes_6Methanotrophs.csv", df_out)
################################################################################################################################################################

################################################################################################################################################################
# 10 Nitrifiers (Ammonia oxidation)

df_ammo = df[(df.ammonia_oxidation.==1) , :]
Genome_size_ammo = convert(Array{Float64}, df_ammo.bp)
V_cell_ammo = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_ammo)
Min_gen_time_ammo = df_ammo.minimum_generation_time_hours
rrn_copies_ammo = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_ammo)
Gram_stain_ammo = repeat(["+"], size(V_cell_ammo,1))

y_DE = [0.097]
y_EM = ones(size(V_cell_ammo,1))
ρ_ps_ammo = zeros(size(V_cell_ammo,1))
r_ms = zeros(size(V_cell_ammo,1))
for i in 1:size(V_cell_ammo,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_ammo[i]], [Min_gen_time_ammo[i]], [Gram_stain_ammo[i]], [rrn_copies_ammo[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_ammo[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_ammo[i]], [Min_gen_time_ammo[i]], [Gram_stain_ammo[i]], [rrn_copies_ammo[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV                           = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_ammo)
r_growth                       = @. (1-1/y_EV)*r_ms
r_main                         = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_ammo, Min_gen_time_ammo, Gram_stain_ammo)
N_SB_ammo                      = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_ammo[i]], ρ_ps_ammo[i].*ones(1,1), [Min_gen_time_ammo[i]], [Gram_stain_ammo[i]])[1] for i in 1:size(V_cell_ammo,1)]
k2p                            = 180.0*60^2
r_assim                        = @. y_DE*N_C*k2p*N_SB_ammo
r_resp_ammo                    = r_growth .+ r_main .+r_assim
median(r_resp_ammo)
std(r_resp_ammo)
# half-saturation constant
Molecular_weight       = [18.04]
D_S                    = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_ammo       = DEBmicroTrait.specific_reference_affinity(V_cell_ammo, reshape(ρ_ps_ammo, 9,1), D_S)
median(K_D_ammo.*14.007)
std(K_D_ammo.*14.007)
K_D_ammo= K_D_ammo.*14.007

df_out = DataFrame()
df_out.genome = vcat(df_ammo.genome)
df_out.domaine = vcat(df_ammo.domain)
df_out.phylum = vcat(df_ammo.phylum)
df_out.habitat_binnedfrom = vcat(df_ammo.habitat_binnedfrom)
df_out.depth_binnedfrom = vcat(df_ammo.depth_binnedfrom)
df_out.functional_group = repeat(["Ammonia oxidizer"], 9)
df_out.substrate = vcat(repeat(["Ammonium"], 9))
df_out.Rmax = vcat(r_resp_ammo)
df_out.K = vcat(K_D_ammo[:,1])
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/emerge/files/TraitPredictionsIsogenieGenomes_10AmmoniaOxidizers.csv", df_out)
################################################################################################################################################################

################################################################################################################################################################
# 11 Nitrifiers (Nitrite oxidation)

df_nitriteox = df[(df.nitrite_oxidation.==1) , :]
Genome_size_nitriteox = convert(Array{Float64}, df_nitriteox.bp)
V_cell_nitriteox = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_nitriteox)
Min_gen_time_nitriteox = df_nitriteox.minimum_generation_time_hours
rrn_copies_nitriteox = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_nitriteox)
Gram_stain_nitriteox = repeat(["+"], size(V_cell_nitriteox,1))

y_DE = [0.0356]
y_EM = ones(size(V_cell_nitriteox,1))
ρ_ps_nitriteox = zeros(size(V_cell_nitriteox,1))
r_ms = zeros(size(V_cell_nitriteox,1))
for i in 1:size(V_cell_nitriteox,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_nitriteox[i]], [Min_gen_time_nitriteox[i]], [Gram_stain_nitriteox[i]], [rrn_copies_nitriteox[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_nitriteox[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_nitriteox[i]], [Min_gen_time_nitriteox[i]], [Gram_stain_nitriteox[i]], [rrn_copies_nitriteox[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV                           = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_nitriteox)
r_growth                       = @. (1-1/y_EV)*r_ms
r_main                         = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_nitriteox, Min_gen_time_nitriteox, Gram_stain_nitriteox)
N_SB_nitriteox                 = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_nitriteox[i]], ρ_ps_nitriteox[i].*ones(1,1), [Min_gen_time_nitriteox[i]], [Gram_stain_nitriteox[i]])[1] for i in 1:size(V_cell_nitriteox,1)]
k2p                            = 180.0*60^2
r_assim                        = @. y_DE*N_C*k2p*N_SB_nitriteox
r_resp_nitriteox               = r_growth .+ r_main .+r_assim
median(r_resp_nitriteox)
std(r_resp_nitriteox)
# half-saturation constant
Molecular_weight       = [46.006]
D_S                    = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_nitriteox       = DEBmicroTrait.specific_reference_affinity(V_cell_nitriteox, reshape(ρ_ps_nitriteox, 33,1), D_S)
median(K_D_nitriteox.*14.007)
std(K_D_nitriteox.*14.007)
K_D_nitriteox= K_D_nitriteox.*14.007

df_out = DataFrame()
df_out.genome = vcat(df_nitriteox.genome)
df_out.domaine = vcat(df_nitriteox.domain)
df_out.phylum = vcat(df_nitriteox.phylum)
df_out.habitat_binnedfrom = vcat(df_nitriteox.habitat_binnedfrom)
df_out.depth_binnedfrom = vcat(df_nitriteox.depth_binnedfrom)
df_out.functional_group = repeat(["Nitrite oxidizer"], 33)
df_out.substrate = vcat(repeat(["Nitrite"], 33))
df_out.Rmax = vcat(r_resp_nitriteox)
df_out.K = vcat(K_D_nitriteox[:,1])
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/emerge/files/TraitPredictionsIsogenieGenomes_11NitriteOxidizers.csv", df_out)
################################################################################################################################################################

################################################################################################################################################################
# 12 Nitrifiers (Nitrite reduction)

df_nitritered = df[(df.nitrite_nitricoxide_partial.==1) , :]
Genome_size_nitritered = convert(Array{Float64}, df_nitritered.bp)
V_cell_nitritered = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_nitritered)
Min_gen_time_nitritered = df_nitritered.minimum_generation_time_hours
rrn_copies_nitritered = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_nitritered)
Gram_stain_nitritered = repeat(["+"], size(V_cell_nitritered,1))

y_DE = [0.118]
y_EM = ones(size(V_cell_nitritered,1))
ρ_ps_nitritered = zeros(size(V_cell_nitritered,1))
r_ms = zeros(size(V_cell_nitritered,1))
for i in 1:size(V_cell_nitritered,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_nitritered[i]], [Min_gen_time_nitritered[i]], [Gram_stain_nitritered[i]], [rrn_copies_nitritered[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_nitritered[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_nitritered[i]], [Min_gen_time_nitritered[i]], [Gram_stain_nitritered[i]], [rrn_copies_nitritered[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV                           = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_nitritered)
r_growth                       = @. (1-1/y_EV)*r_ms
r_main                         = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_nitritered, Min_gen_time_nitritered, Gram_stain_nitritered)
N_SB_nitritered                 = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_nitritered[i]], ρ_ps_nitritered[i].*ones(1,1), [Min_gen_time_nitritered[i]], [Gram_stain_nitritered[i]])[1] for i in 1:size(V_cell_nitritered,1)]
k2p                            = 180.0*60^2
r_assim                        = @. y_DE*N_C*k2p*N_SB_nitritered
r_resp_nitritered               = r_growth .+ r_main .+r_assim
median(r_resp_nitritered)
std(r_resp_nitritered)
# half-saturation constant
Molecular_weight       = [18.04]
D_S                    = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_nitritered       = DEBmicroTrait.specific_reference_affinity(V_cell_nitritered, reshape(ρ_ps_nitritered, 40,1), D_S)
median(K_D_nitritered.*14.007)
std(K_D_nitritered.*14.007)
K_D_nitritered= K_D_nitritered.*14.007

df_out = DataFrame()
df_out.genome = vcat(df_nitritered.genome)
df_out.domaine = vcat(df_nitritered.domain)
df_out.phylum = vcat(df_nitritered.phylum)
df_out.habitat_binnedfrom = vcat(df_nitritered.habitat_binnedfrom)
df_out.depth_binnedfrom = vcat(df_nitritered.depth_binnedfrom)
df_out.functional_group = repeat(["Nitrite reducers"], 40)
df_out.substrate = vcat(repeat(["Nitrite"], 40))
df_out.Rmax = vcat(r_resp_nitritered)
df_out.K = vcat(K_D_nitritered[:,1])
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/emerge/files/TraitPredictionsIsogenieGenomes_12NitriteReducers.csv", df_out)
################################################################################################################################################################

################################################################################################################################################################
# 7 Denitrifiers (Nitrate reductase)

df_nitratered = df[(df.nitrite_nitrate_partial.==1) , :]
Genome_size_nitratered = convert(Array{Float64}, df_nitratered.bp)
V_cell_nitratered = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_nitratered)
Min_gen_time_nitratered = df_nitratered.minimum_generation_time_hours
rrn_copies_nitratered = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_nitratered)
Gram_stain_nitratered = repeat(["+"], size(V_cell_nitratered,1))

N_C  = [6]
y_DE = [2.4]./N_C
y_EM = ones(size(V_cell_nitratered,1))
ρ_ps_nitratered = zeros(size(V_cell_nitratered,1))
r_ms = zeros(size(V_cell_nitratered,1))
for i in 1:size(V_cell_nitratered,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_nitratered[i]], [Min_gen_time_nitratered[i]], [Gram_stain_nitratered[i]], [rrn_copies_nitratered[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_nitratered[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_nitratered[i]], [Min_gen_time_nitratered[i]], [Gram_stain_nitratered[i]], [rrn_copies_nitratered[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV                           = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_nitratered)
r_growth                       = @. (1-1/y_EV)*r_ms
r_main                         = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_nitratered, Min_gen_time_nitratered, Gram_stain_nitratered)
N_SB_nitratered                 = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_nitratered[i]], ρ_ps_nitratered[i].*ones(1,1), [Min_gen_time_nitratered[i]], [Gram_stain_nitratered[i]])[1] for i in 1:size(V_cell_nitratered,1)]
k2p                            = 180.0*60^2
r_assim                        = @. y_DE*N_C*k2p*N_SB_nitratered
r_resp_nitratered               = r_growth .+ r_main .+r_assim
median(r_resp_nitratered)
std(r_resp_nitratered)
# half-saturation constant
Molecular_weight       = [180.16]
D_S                    = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_nitratered       = DEBmicroTrait.specific_reference_affinity(V_cell_nitratered, reshape(ρ_ps_nitratered, 40,1), D_S)
median(K_D_nitratered.*N_C*12.011)
std(K_D_nitratered.*N_C*12.011)
K_D_nitratered= K_D_nitratered.*N_C*12.011

df_out = DataFrame()
df_out.genome = vcat(df_nitratered.genome)
df_out.domaine = vcat(df_nitratered.domain)
df_out.phylum = vcat(df_nitratered.phylum)
df_out.habitat_binnedfrom = vcat(df_nitratered.habitat_binnedfrom)
df_out.depth_binnedfrom = vcat(df_nitratered.depth_binnedfrom)
df_out.functional_group = repeat(["Nitrate reductase"], 40)
df_out.substrate = vcat(repeat(["Nitrate"], 40))
df_out.Rmax = vcat(r_resp_nitratered)
df_out.K = vcat(K_D_nitratered[:,1])
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/emerge/files/TraitPredictionsIsogenieGenomes_7NitrateReductase.csv", df_out)
################################################################################################################################################################

################################################################################################################################################################
# 8 Denitrifiers (Nitrite reductase)

df_nitritereductase = df[(df.nitrite_nitricoxide_partial.==1) , :]
Genome_size_nitritereductase = convert(Array{Float64}, df_nitritereductase.bp)
V_cell_nitritereductase = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_nitritereductase)
Min_gen_time_nitritereductase = df_nitritereductase.minimum_generation_time_hours
rrn_copies_nitritereductase = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_nitritereductase)
Gram_stain_nitritereductase = repeat(["+"], size(V_cell_nitritereductase,1))

N_C  = [6]
y_DE = [3.11]./N_C
y_EM = ones(size(V_cell_nitritereductase,1))
ρ_ps_nitritereductase = zeros(size(V_cell_nitritereductase,1))
r_ms = zeros(size(V_cell_nitritereductase,1))
for i in 1:size(V_cell_nitritereductase,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_nitritereductase[i]], [Min_gen_time_nitritereductase[i]], [Gram_stain_nitritereductase[i]], [rrn_copies_nitritereductase[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_nitritereductase[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_nitritereductase[i]], [Min_gen_time_nitritereductase[i]], [Gram_stain_nitritereductase[i]], [rrn_copies_nitritereductase[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV                           = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_nitritereductase)
r_growth                       = @. (1-1/y_EV)*r_ms
r_main                         = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_nitritereductase, Min_gen_time_nitritereductase, Gram_stain_nitritereductase)
N_SB_nitritereductase                 = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_nitritereductase[i]], ρ_ps_nitritereductase[i].*ones(1,1), [Min_gen_time_nitritereductase[i]], [Gram_stain_nitritereductase[i]])[1] for i in 1:size(V_cell_nitritereductase,1)]
k2p                            = 180.0*60^2
r_assim                        = @. y_DE*N_C*k2p*N_SB_nitritereductase
r_resp_nitritereductase               = r_growth .+ r_main .+r_assim
median(r_resp_nitritereductase)
std(r_resp_nitritereductase)
# half-saturation constant
Molecular_weight       = [180.16]
D_S                    = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_nitritereductase       = DEBmicroTrait.specific_reference_affinity(V_cell_nitritereductase, reshape(ρ_ps_nitritereductase, 40,1), D_S)
median(K_D_nitritereductase.*N_C*12.011)
std(K_D_nitritereductase.*N_C*12.011)
K_D_nitritereductase= K_D_nitritereductase.*N_C*12.011

df_out = DataFrame()
df_out.genome = vcat(df_nitritereductase.genome)
df_out.domaine = vcat(df_nitritereductase.domain)
df_out.phylum = vcat(df_nitritereductase.phylum)
df_out.habitat_binnedfrom = vcat(df_nitritereductase.habitat_binnedfrom)
df_out.depth_binnedfrom = vcat(df_nitritereductase.depth_binnedfrom)
df_out.functional_group = repeat(["Nitrite reductase"], 40)
df_out.substrate = vcat(repeat(["Nitrite"], 40))
df_out.Rmax = vcat(r_resp_nitritereductase)
df_out.K = vcat(K_D_nitritereductase[:,1])
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/emerge/files/TraitPredictionsIsogenieGenomes_8NitriteReductase.csv", df_out)

################################################################################################################################################################

################################################################################################################################################################
# 9 Denitrifiers (Nitrous oxide reductase)

df_nitrousoxidereductase = df[(df.nitrousoxide_dinitrogen.==1) , :]
Genome_size_nitrousoxidereductase = convert(Array{Float64}, df_nitrousoxidereductase.bp)
V_cell_nitrousoxidereductase = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_nitrousoxidereductase)
Min_gen_time_nitrousoxidereductase = df_nitrousoxidereductase.minimum_generation_time_hours
rrn_copies_nitrousoxidereductase = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_nitrousoxidereductase)
Gram_stain_nitrousoxidereductase = repeat(["+"], size(V_cell_nitrousoxidereductase,1))

N_C  = [6]
y_DE = [3.52]./N_C
y_EM = ones(size(V_cell_nitrousoxidereductase,1))
ρ_ps_nitrousoxidereductase = zeros(size(V_cell_nitrousoxidereductase,1))
r_ms = zeros(size(V_cell_nitrousoxidereductase,1))
for i in 1:size(V_cell_nitrousoxidereductase,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_nitrousoxidereductase[i]], [Min_gen_time_nitrousoxidereductase[i]], [Gram_stain_nitrousoxidereductase[i]], [rrn_copies_nitrousoxidereductase[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_nitrousoxidereductase[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_nitrousoxidereductase[i]], [Min_gen_time_nitrousoxidereductase[i]], [Gram_stain_nitrousoxidereductase[i]], [rrn_copies_nitrousoxidereductase[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV                           = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_nitrousoxidereductase)
r_growth                       = @. (1-1/y_EV)*r_ms
r_main                         = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_nitrousoxidereductase, Min_gen_time_nitrousoxidereductase, Gram_stain_nitrousoxidereductase)
N_SB_nitrousoxidereductase                 = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_nitrousoxidereductase[i]], ρ_ps_nitrousoxidereductase[i].*ones(1,1), [Min_gen_time_nitrousoxidereductase[i]], [Gram_stain_nitrousoxidereductase[i]])[1] for i in 1:size(V_cell_nitrousoxidereductase,1)]
k2p                            = 180.0*60^2
r_assim                        = @. y_DE*N_C*k2p*N_SB_nitrousoxidereductase
r_resp_nitrousoxidereductase               = r_growth .+ r_main .+r_assim
median(r_resp_nitrousoxidereductase)
std(r_resp_nitrousoxidereductase)
# half-saturation constant
Molecular_weight       = [180.16]
D_S                    = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_nitrousoxidereductase       = DEBmicroTrait.specific_reference_affinity(V_cell_nitrousoxidereductase, reshape(ρ_ps_nitrousoxidereductase, 16,1), D_S)
median(K_D_nitrousoxidereductase.*N_C*12.011)
std(K_D_nitrousoxidereductase.*N_C*12.011)
K_D_nitrousoxidereductase= K_D_nitrousoxidereductase.*N_C*12.011

df_out = DataFrame()
df_out.genome = vcat(df_nitrousoxidereductase.genome)
df_out.domaine = vcat(df_nitrousoxidereductase.domain)
df_out.phylum = vcat(df_nitrousoxidereductase.phylum)
df_out.habitat_binnedfrom = vcat(df_nitrousoxidereductase.habitat_binnedfrom)
df_out.depth_binnedfrom = vcat(df_nitrousoxidereductase.depth_binnedfrom)
df_out.functional_group = repeat(["Nitrous oxide reductase"], 16)
df_out.substrate = vcat(repeat(["Nitrous oxide"], 16))
df_out.Rmax = vcat(r_resp_nitrousoxidereductase)
df_out.K = vcat(K_D_nitrousoxidereductase[:,1])
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/emerge/files/TraitPredictionsIsogenieGenomes_9NitrousoxideReductase.csv", df_out)

################################################################################################################################################################

################################################################################################################################################################
# 4 Acetoclastic methanogens

df_acetoclastic = df[(df.acetoclastic_methanogenesis.==1) , :]
Genome_size_acetoclastic = convert(Array{Float64}, df_acetoclastic.bp)
V_cell_acetoclastic = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_acetoclastic)
Min_gen_time_acetoclastic = df_acetoclastic.minimum_generation_time_hours
rrn_copies_acetoclastic = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_acetoclastic)
Gram_stain_acetoclastic = repeat(["+"], size(V_cell_acetoclastic,1))

N_C  = [2]
y_DE = [0.06]./N_C
y_EM = ones(size(V_cell_acetoclastic,1))
ρ_ps_acetoclastic = zeros(size(V_cell_acetoclastic,1))
r_ms = zeros(size(V_cell_acetoclastic,1))
for i in 1:size(V_cell_acetoclastic,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_acetoclastic[i]], [Min_gen_time_acetoclastic[i]], [Gram_stain_acetoclastic[i]], [rrn_copies_acetoclastic[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_acetoclastic[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_acetoclastic[i]], [Min_gen_time_acetoclastic[i]], [Gram_stain_acetoclastic[i]], [rrn_copies_acetoclastic[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV                           = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_acetoclastic)
r_growth                       = @. (1-1/y_EV)*r_ms
r_main                         = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_acetoclastic, Min_gen_time_acetoclastic, Gram_stain_acetoclastic)
N_SB_acetoclastic              = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_acetoclastic[i]], ρ_ps_acetoclastic[i].*ones(1,1), [Min_gen_time_acetoclastic[i]], [Gram_stain_acetoclastic[i]])[1] for i in 1:size(V_cell_acetoclastic,1)]
k2p                            = 180.0*60^2
r_assim                        = @. y_DE*N_C*k2p*N_SB_acetoclastic
r_resp_acetoclastic               = r_growth .+ r_main .+r_assim
median(r_resp_acetoclastic)
std(r_resp_acetoclastic)
# half-saturation constant
Molecular_weight       = [60.05]
D_S                    = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_acetoclastic       = DEBmicroTrait.specific_reference_affinity(V_cell_acetoclastic, reshape(ρ_ps_acetoclastic*100, 17,1), D_S)
median(K_D_acetoclastic.*N_C*12.011)
std(K_D_acetoclastic.*N_C*12.011)
K_D_acetoclastic= K_D_acetoclastic.*N_C*12.011

df_out = DataFrame()
df_out.genome = vcat(df_acetoclastic.genome)
df_out.domaine = vcat(df_acetoclastic.domain)
df_out.phylum = vcat(df_acetoclastic.phylum)
df_out.habitat_binnedfrom = vcat(df_acetoclastic.habitat_binnedfrom)
df_out.depth_binnedfrom = vcat(df_acetoclastic.depth_binnedfrom)
df_out.functional_group = repeat(["Acetoclastic methanogen"], 17)
df_out.substrate = vcat(repeat(["Acetic acid"], 17))
df_out.Rmax = vcat(r_resp_acetoclastic)
df_out.K = vcat(K_D_acetoclastic[:,1])
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/emerge/files/TraitPredictionsIsogenieGenomes_4AcetoclasticMethanogens.csv", df_out)

################################################################################################################################################################

################################################################################################################################################################
# 5 Hydrogenotrophic methanogens

df_hydrogenotrophic = df[(df.hydrogenotrophic_methanogenesis.==1) , :]
Genome_size_hydrogenotrophic = convert(Array{Float64}, df_hydrogenotrophic.bp)
V_cell_hydrogenotrophic = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_hydrogenotrophic)
Min_gen_time_hydrogenotrophic = df_hydrogenotrophic.minimum_generation_time_hours
rrn_copies_hydrogenotrophic = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_hydrogenotrophic)
Gram_stain_hydrogenotrophic = repeat(["+"], size(V_cell_hydrogenotrophic,1))

N_C  = [1]
y_DE = [0.24]./N_C
y_EM = ones(size(V_cell_hydrogenotrophic,1))
ρ_ps_hydrogenotrophic = zeros(size(V_cell_hydrogenotrophic,1))
r_ms = zeros(size(V_cell_hydrogenotrophic,1))
for i in 1:size(V_cell_hydrogenotrophic,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_hydrogenotrophic[i]], [Min_gen_time_hydrogenotrophic[i]], [Gram_stain_hydrogenotrophic[i]], [rrn_copies_hydrogenotrophic[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_hydrogenotrophic[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_hydrogenotrophic[i]], [Min_gen_time_hydrogenotrophic[i]], [Gram_stain_hydrogenotrophic[i]], [rrn_copies_hydrogenotrophic[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV                           = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_hydrogenotrophic)
r_growth                       = @. (1-1/y_EV)*r_ms
r_main                         = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_hydrogenotrophic, Min_gen_time_hydrogenotrophic, Gram_stain_hydrogenotrophic)
N_SB_hydrogenotrophic              = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_hydrogenotrophic[i]], ρ_ps_hydrogenotrophic[i].*ones(1,1), [Min_gen_time_hydrogenotrophic[i]], [Gram_stain_hydrogenotrophic[i]])[1] for i in 1:size(V_cell_hydrogenotrophic,1)]
k2p                            = 180.0*60^2
r_assim                        = @. y_DE*N_C*k2p*N_SB_hydrogenotrophic
r_resp_hydrogenotrophic               = r_growth .+ r_main .+r_assim
median(r_resp_hydrogenotrophic)
std(r_resp_hydrogenotrophic)
# half-saturation constant
Molecular_weight       = [44.009]
D_S                    = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_hydrogenotrophic       = DEBmicroTrait.specific_reference_affinity(V_cell_hydrogenotrophic, reshape(ρ_ps_hydrogenotrophic, 18,1), D_S)
median(K_D_hydrogenotrophic.*N_C*12.011)
std(K_D_hydrogenotrophic.*N_C*12.011)
K_D_hydrogenotrophic= K_D_hydrogenotrophic.*N_C*12.011

df_out = DataFrame()
df_out.genome = vcat(df_hydrogenotrophic.genome)
df_out.domaine = vcat(df_hydrogenotrophic.domain)
df_out.phylum = vcat(df_hydrogenotrophic.phylum)
df_out.habitat_binnedfrom = vcat(df_hydrogenotrophic.habitat_binnedfrom)
df_out.depth_binnedfrom = vcat(df_hydrogenotrophic.depth_binnedfrom)
df_out.functional_group = repeat(["Hydrogenotrophic methanogen"], 18)
df_out.substrate = vcat(repeat(["Carbon dioxide"], 18))
df_out.Rmax = vcat(r_resp_hydrogenotrophic)
df_out.K = vcat(K_D_hydrogenotrophic[:,1])
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/emerge/files/TraitPredictionsIsogenieGenomes_5HydrogenotrophicMethanogens.csv", df_out)

################################################################################################################################################################
# 13 Aerobic diazotroph

df_aerobicdiaz = df[(df.aerobic_diazotroph.==1) , :]
Genome_size_aerobicdiaz = convert(Array{Float64}, df_aerobicdiaz.bp)
V_cell_aerobicdiaz = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_aerobicdiaz)
Min_gen_time_aerobicdiaz = df_aerobicdiaz.minimum_generation_time_hours
rrn_copies_aerobicdiaz = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_aerobicdiaz)
Gram_stain_aerobicdiaz = repeat(["+"], size(V_cell_aerobicdiaz,1))

# Glucose
N_C  = [6]
y_DE = [1.29]./N_C
y_EM = ones(size(V_cell_aerobicdiaz,1))
ρ_ps_aerobicdiaz = zeros(size(V_cell_aerobicdiaz,1))
r_ms = zeros(size(V_cell_aerobicdiaz,1))
for i in 1:size(V_cell_aerobicdiaz,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_aerobicdiaz[i]], [Min_gen_time_aerobicdiaz[i]], [Gram_stain_aerobicdiaz[i]], [rrn_copies_aerobicdiaz[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_aerobicdiaz[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_aerobicdiaz[i]], [Min_gen_time_aerobicdiaz[i]], [Gram_stain_aerobicdiaz[i]], [rrn_copies_aerobicdiaz[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV                           = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_aerobicdiaz)
r_growth                       = @. (1-1/y_EV)*r_ms
r_main                         = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_aerobicdiaz, Min_gen_time_aerobicdiaz, Gram_stain_aerobicdiaz)
N_SB_aerobicdiaz              = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_aerobicdiaz[i]], ρ_ps_aerobicdiaz[i].*ones(1,1), [Min_gen_time_aerobicdiaz[i]], [Gram_stain_aerobicdiaz[i]])[1] for i in 1:size(V_cell_aerobicdiaz,1)]
k2p                            = 180.0*60^2
r_assim                        = @. y_DE*N_C*k2p*N_SB_aerobicdiaz
r_resp_aerobicdiaz               = r_growth .+ r_main .+r_assim
median(r_resp_aerobicdiaz)
std(r_resp_aerobicdiaz)
# half-saturation constant
Molecular_weight       = [180.16]
D_S                    = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_aerobicdiaz       = DEBmicroTrait.specific_reference_affinity(V_cell_aerobicdiaz, reshape(ρ_ps_aerobicdiaz, 20,1), D_S)
median(K_D_aerobicdiaz.*N_C*12.011)
std(K_D_aerobicdiaz.*N_C*12.011)
K_D_aerobicdiaz= K_D_aerobicdiaz.*N_C*12.011

df_out = DataFrame()
df_out.genome = vcat(df_aerobicdiaz.genome)
df_out.domaine = vcat(df_aerobicdiaz.domain)
df_out.phylum = vcat(df_aerobicdiaz.phylum)
df_out.habitat_binnedfrom = vcat(df_aerobicdiaz.habitat_binnedfrom)
df_out.depth_binnedfrom = vcat(df_aerobicdiaz.depth_binnedfrom)
df_out.functional_group = repeat(["Aerobic diazotroph"], 20)
df_out.substrate = vcat(repeat(["Glucose"], 20))
df_out.Rmax = vcat(r_resp_aerobicdiaz)
df_out.K = vcat(K_D_aerobicdiaz[:,1])
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/emerge/files/TraitPredictionsIsogenieGenomes_13AerobicDiazotroph.csv", df_out)

################################################################################################################################################################
# 14 Anaerobic diazotroph

df_anaerobicdiaz = df[(df.aerobic_diazotroph.==1) , :]
Genome_size_anaerobicdiaz = convert(Array{Float64}, df_anaerobicdiaz.bp)
V_cell_anaerobicdiaz = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_anaerobicdiaz)
Min_gen_time_anaerobicdiaz = df_anaerobicdiaz.minimum_generation_time_hours
rrn_copies_anaerobicdiaz = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_anaerobicdiaz)
Gram_stain_anaerobicdiaz = repeat(["+"], size(V_cell_anaerobicdiaz,1))

# CO2
N_C  = [1]
y_DE = [0.27]./N_C
y_EM = ones(size(V_cell_anaerobicdiaz,1))
ρ_ps_anaerobicdiaz = zeros(size(V_cell_anaerobicdiaz,1))
r_ms = zeros(size(V_cell_anaerobicdiaz,1))
for i in 1:size(V_cell_anaerobicdiaz,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_anaerobicdiaz[i]], [Min_gen_time_anaerobicdiaz[i]], [Gram_stain_anaerobicdiaz[i]], [rrn_copies_anaerobicdiaz[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_anaerobicdiaz[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_anaerobicdiaz[i]], [Min_gen_time_anaerobicdiaz[i]], [Gram_stain_anaerobicdiaz[i]], [rrn_copies_anaerobicdiaz[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV                           = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_anaerobicdiaz)
r_growth                       = @. (1-1/y_EV)*r_ms
r_main                         = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_anaerobicdiaz, Min_gen_time_anaerobicdiaz, Gram_stain_anaerobicdiaz)
N_SB_anaerobicdiaz              = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_anaerobicdiaz[i]], ρ_ps_anaerobicdiaz[i].*ones(1,1), [Min_gen_time_anaerobicdiaz[i]], [Gram_stain_anaerobicdiaz[i]])[1] for i in 1:size(V_cell_anaerobicdiaz,1)]
k2p                            = 180.0*60^2
r_assim                        = @. y_DE*N_C*k2p*N_SB_anaerobicdiaz
r_resp_anaerobicdiaz               = r_growth .+ r_main .+r_assim
median(r_resp_anaerobicdiaz)
std(r_resp_anaerobicdiaz)
# half-saturation constant
Molecular_weight       = [44.009]
D_S                    = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_anaerobicdiaz       = DEBmicroTrait.specific_reference_affinity(V_cell_anaerobicdiaz, reshape(ρ_ps_anaerobicdiaz, 20,1), D_S)
median(K_D_anaerobicdiaz.*N_C*12.011)
std(K_D_anaerobicdiaz.*N_C*12.011)
K_D_anaerobicdiaz= K_D_anaerobicdiaz.*N_C*12.011

df_out = DataFrame()
df_out.genome = vcat(df_anaerobicdiaz.genome)
df_out.domaine = vcat(df_anaerobicdiaz.domain)
df_out.phylum = vcat(df_anaerobicdiaz.phylum)
df_out.habitat_binnedfrom = vcat(df_anaerobicdiaz.habitat_binnedfrom)
df_out.depth_binnedfrom = vcat(df_anaerobicdiaz.depth_binnedfrom)
df_out.functional_group = repeat(["Anaerobic diazotroph"], 20)
df_out.substrate = vcat(repeat(["Carbon dioxide"], 20))
df_out.Rmax = vcat(r_resp_anaerobicdiaz)
df_out.K = vcat(K_D_anaerobicdiaz[:,1])
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/emerge/files/TraitPredictionsIsogenieGenomes_14AnaerobicDiazotroph.csv", df_out)
