using UnROOT, SHA, Printf

"""
    prepare_cms(path)

Read all 61,540,413 events from CERN Open Data record 12341 (CC0-1.0):
Wunsch (2019), DOI 10.7483/OPENDATA.CMS.LVG5.QT81.
Input: Run2012BC_DoubleMuParked_Muons.root, 2,244,449,133 bytes.
Requires UnROOT.jl; neither ROOT nor ScientificFitting is needed.

Select exactly two oppositely charged muons, as in ROOT's df102 tutorial.
Return every reconstructed mass in [2.8, 3.4) GeV/c^2 for unbinned fitting,
the two muon transverse momenta and pseudorapidities, pair transverse momentum
and rapidity, cut counts, and a logarithmically binned overview spectrum.
Kinematic arrays have one entry per retained mass, in the same order. No sampling.
The SHA-256 was obtained after verifying CERN's Adler-32 checksum 1fa61aca.
"""
function prepare_cms(path)
    @assert filesize(path) == 2244449133 "Wrong or incomplete ROOT file"
    @assert bytes2hex(open(sha256, path)) ==
        "f8a9d40dba9ee7131a4110f93b4e14a29f99255cbff777018848477f2f61b00b"
    events = LazyTree(ROOTFile(path), "Events",
        ["nMuon", "Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_charge"])
    @assert length(events) == 61540413
    two = opposite = invalid = 0
    mass = Float64[]
    muon_pt, muon_eta = NTuple{2,Float64}[], NTuple{2,Float64}[]
    pair_pt, rapidity = Float64[], Float64[]
    for values in (mass, muon_pt, muon_eta, pair_pt, rapidity)
        sizehint!(values, 3_000_000)
    end
    edges = exp10.(range(log10(0.25), log10(300.); length=301))
    counts = zeros(Int, 300)
    for event in events
        event.nMuon == 2 || continue
        two += 1
        event.Muon_charge[1] != event.Muon_charge[2] || continue
        opposite += 1
        # Promote Float32 inputs before four-vector arithmetic, using c=1.
        pt = ntuple(j -> Float64(event.Muon_pt[j]), 2)
        eta = ntuple(j -> Float64(event.Muon_eta[j]), 2)
        phi = ntuple(j -> Float64(event.Muon_phi[j]), 2)
        m = ntuple(j -> Float64(event.Muon_mass[j]), 2)
        pz = pt .* sinh.(eta)
        energy = sqrt.(pt.^2 .+ pz.^2 .+ m.^2)
        m2 = m[1]^2 + m[2]^2 + 2*(energy[1]*energy[2] - pz[1]*pz[2] -
            pt[1]*pt[2]*cos(phi[1]-phi[2]))
        if !(isfinite(m2) && m2 >= 0)
            invalid += 1
            continue
        end
        value = sqrt(m2)
        bin = searchsortedlast(edges, value)
        1 <= bin <= length(counts) && (counts[bin] += 1)
        if 2.8 <= value < 3.4
            # Retain the same candidates for resolution checks versus kinematics.
            pt_pair = hypot(pt[1]*cos(phi[1]) + pt[2]*cos(phi[2]),
                pt[1]*sin(phi[1]) + pt[2]*sin(phi[2]))
            push!(mass, value)
            push!(muon_pt, pt)
            push!(muon_eta, eta)
            push!(pair_pt, pt_pair)
            push!(rapidity, asinh(sum(pz)/hypot(value, pt_pair)))
        end
    end
    @assert two == 31104343 && opposite == 24067843 && invalid == 0
    @assert length(mass) == 2551454
    return (; events=length(events), two, opposite, invalid, mass, muon_pt,
        muon_eta, pair_pt, rapidity, edges, counts)
end

"""
    write_cms_histogram(path, data)

Export 300 equal mass bins and three disjoint max(abs(muon_eta)) groups:
[0, 0.9), [0.9, 1.4), [1.4, Inf). No candidates are discarded by grouping.
The CSV is the compact, reproducible input used by the documentation build.
"""
function write_cms_histogram(path, data)
    edges = collect(range(2.8, 3.4; length=301))
    eta_edges = [0., 0.9, 1.4, Inf]
    counts = zeros(Int, 300, 3)
    @assert length(data.mass) == length(data.muon_eta)
    for (mass, eta) in zip(data.mass, data.muon_eta)
        counts[searchsortedlast(edges, mass),
            searchsortedlast(eta_edges, maximum(abs, eta))] += 1
    end
    @assert sum(counts) == length(data.mass) == 2_551_454
    open(path, "w") do io
        println(io, "left_GeV,right_GeV,total,eta_0_09,eta_09_14,eta_14_inf")
        for i in axes(counts, 1)
            @printf(io, "%.3f,%.3f,%d,%d,%d,%d\n", edges[i], edges[i+1],
                sum(counts[i,:]), counts[i,1], counts[i,2], counts[i,3])
        end
    end
    return path
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("Usage: julia prepare.jl Run2012BC_DoubleMuParked_Muons.root")
    data = prepare_cms(only(ARGS))
    println("Events: ", data.events, "; exactly two muons: ", data.two,
        "; opposite charge: ", data.opposite, "; J/psi window: ", length(data.mass))
end
