"""
    prepare_lhcb(path)

Rebuild the unweighted KKK mass histogram used in the LHCb gallery page.
Input: CERN Open Data record 4900, B2HHH_MagnetUp.root (CC0-1.0).
The SHA-256 below was computed after verifying CERN's Adler-32 0568274c.
Requires UnROOT.jl 0.11.11; neither ROOT nor the fitting package is needed.
"""
function prepare_lhcb(path)
    @assert filesize(path) == 444723234 "Wrong or incomplete ROOT file"
    @assert bytes2hex(open(sha256, path)) ==
        "c42ad9e47931e1404bf94ad82ea22a0acd10bc9cfbb58e77a6b0fff08ead7859"
    tree = LazyTree(ROOTFile(path), "DecayTree",
        [r"H[123]_(PX|PY|PZ|ProbK|ProbPi|isMuon|Charge)$"])
    edges, counts = collect(5200.:5.:5600.), zeros(Int, 80)
    selected = invalid = below = above = 0
    for event in tree
        # Use the official notebook's particle-ID cuts, without a charge split.
        event.H1_ProbK > 0.5 && event.H2_ProbK > 0.5 && event.H3_ProbK > 0.5 || continue
        event.H1_ProbPi < 0.5 && event.H2_ProbPi < 0.5 && event.H3_ProbPi < 0.5 || continue
        iszero(event.H1_isMuon) && iszero(event.H2_isMuon) && iszero(event.H3_isMuon) || continue
        selected += 1
        @assert abs(event.H1_Charge + event.H2_Charge + event.H3_Charge) == 1
        px = (event.H1_PX, event.H2_PX, event.H3_PX)
        py = (event.H1_PY, event.H2_PY, event.H3_PY)
        pz = (event.H1_PZ, event.H2_PZ, event.H3_PZ)
        # Momenta are in MeV/c. Assign the charged-kaon mass to each track;
        # use c=1 while summing four-momenta, then report m in MeV/c^2.
        energy = sum(sqrt(px[i]^2 + py[i]^2 + pz[i]^2 + 493.677^2) for i in 1:3)
        m2 = energy^2 - sum(px)^2 - sum(py)^2 - sum(pz)^2
        if !isfinite(m2) || m2 <= 0
            invalid += 1
            continue
        end
        mass = sqrt(m2)
        if mass < first(edges)
            below += 1
        elseif mass >= last(edges)
            above += 1
        else
            counts[searchsortedlast(edges, mass)] += 1 # [lower, upper) bins
        end
    end
    @assert selected == invalid + below + above + sum(counts)
    @assert invalid == 0
    return (; candidates=length(tree), selected, below, above, edges, counts)
end

using UnROOT, SHA
if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("Usage: julia prepare.jl B2HHH_MagnetUp.root")
    result = prepare_lhcb(only(ARGS))
    println("Input candidates: ", result.candidates, "; selected: ", result.selected)
    println("Below/above fit window: ", (result.below, result.above))
    println("Histogram counts (sum = ", sum(result.counts), "):")
    show(stdout, result.counts)
    println()
end
