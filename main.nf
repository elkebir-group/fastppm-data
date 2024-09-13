params.simulate = 'simulate.py'

params.nmutations = [100, 500, 1000, 2500]
params.nsamples   = [3]
params.seeds      = 1..20
params.coverage   = [30, 100, 1000]

process create_sim {
    cpus 1
    memory '4 GB'
    time '59m'
    publishDir "simulations/n${mutations}_s${samples}_c${coverage}_r${seed}", mode: 'copy', overwrite: true

    input:
        tuple val(mutations), val(samples), val(coverage), val(seed)

    output:
        tuple file("sim_clonal_matrix.txt"), file("sim_frequency_matrix.txt"), 
              file("sim_total_matrix.txt"), file("sim_tree.txt"), file("sim_usage_matrix.txt"), 
              file("sim_variant_matrix.txt"),
              val("n${mutations}_s${samples}_c${coverage}_r${seed}")

    """
    '${params.simulate}' --mutations ${mutations} --samples ${samples} --coverage ${coverage} --seed ${seed} --output sim
    """
}

workflow {
    parameter_channel = channel.fromList(params.nmutations)
                               .combine(channel.fromList(params.nsamples))
                               .combine(channel.fromList(params.coverage))
                               .combine(channel.fromList(params.seeds))

    simulation = parameter_channel | create_sim
}
