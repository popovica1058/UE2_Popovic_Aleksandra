
# Benoetigte Pakete einbinden
using Random

### Datentypen

# Datentypen definieren
struct Instance
    name::String
    n::Int
    x::Vector{Int}
    y::Vector{Int}
    dist::Matrix{Int}
    bestobj::Int

    function Instance(name::String, x::Vector{Int}, y::Vector{Int}, dist::Matrix{Int}, bestopt::Int)
        if length(x) != length(y) || length(x) <= 0
            throw(ArgumentError("x und y muessen gleich lang und nicht leer sein."))
        end
        if size(dist) != (length(x), length(y))
            throw(ArgumentError("dist muss eine n x n Matrix sein."))
        end
        new(name, length(x), x, y, dist, bestopt)
    end
end

mutable struct Solution
    const inst::Instance
    tour::Vector{Int}
    obj::Int

    function Solution(inst::Instance)
        new(inst, Int[], 0)
    end
end

### Beispiel 1

function distmatrix(x::Vector{Int}, y::Vector{Int}; p::Real = 2)::Matrix{Int}
    if p <= 0
        throw(DomainError(p, "p > 0 erwartet!"))
    end
    n = length(x)
    M = zeros(Int, n, n)
    for i in 1:n
        for j in i+1:n
            dist = ((abs(x[i] - x[j])^p + abs(y[i] - y[j])^p))^(1/p)
            M[i, j] = M[j, i] = round(Int, dist)
        end
    end
    return M
end

### Beispiel 2

function nearest_neighbor(inst::Instance, startcity::Int = rand(1:inst.n))::Solution
    n = inst.n
    tour = [startcity]
    unvisited = setdiff(1:n, tour)
    current = startcity

    while !isempty(unvisited)
        next_city = argmin(inst.dist[current, u] for u in unvisited)
        city = unvisited[next_city]
        push!(tour, city)
        deleteat!(unvisited, next_city)
        current = city
    end

    obj = sum(inst.dist[tour[i], tour[i+1]] for i in 1:n-1) + inst.dist[tour[end], tour[1]]
    return Solution(inst, tour, obj)
end

### Beispiel 3

# Hilfsfunktion fuer Cheapest Insertion
function diff_insertion(sol::Solution, pos::Int, city::Int)::Int
    tour = sol.tour
    inst = sol.inst
    n = length(tour)
    if n == 0
        return 0
    elseif n == 1
        return 2 * inst.dist[tour[1], city]
    else
        i = pos == 1 ? n : pos - 1
        j = pos > n ? 1 : pos
        return inst.dist[tour[i], city] + inst.dist[city, tour[j]] - inst.dist[tour[i], tour[j]]
    end
end

# Cheapest Insertion Heuristik
function cheapest_insertion(inst::Instance, order::Vector{Int} = shuffle(1:inst.n))::Solution
    if sort(order) != collect(1:inst.n)
        throw(ArgumentError("order muss eine Permutation von 1:n sein."))
    end
    sol = Solution(inst)
    push!(sol.tour, order[1])
    push!(sol.tour, order[2])
    sol.obj = inst.dist[order[1], order[2]] + inst.dist[order[2], order[1]]

    for city in order[3:end]
        costs = [diff_insertion(sol, pos, city) for pos in 1:length(sol.tour)+1]
        bestpos = argmin(costs)
        insert!(sol.tour, bestpos, city)
        sol.obj += costs[bestpos]
    end

    return sol
end

### Beispiel 4

# Hilfsfunktion fuer 2-Opt
function diff_twoopt(sol::Solution, a::Int, b::Int)::Int
    n = length(sol.tour)
    if a == b
        return 0
    end
    a1 = mod1(a - 1, n)
    b1 = mod1(b + 1, n)
    tour = sol.tour
    inst = sol.inst

    delta = -inst.dist[tour[a1], tour[a]] - inst.dist[tour[b], tour[b1]] +
            inst.dist[tour[a1], tour[b]] + inst.dist[tour[a], tour[b1]]
    return delta
end

# Lokale Suche mit 2-Opt
function localsearch!(sol::Solution)::Nothing
    improved = true
    while improved
        improved = false
        for i in 2:(length(sol.tour)-2)
            for j in (i+1):(length(sol.tour)-1)
                delta = diff_twoopt(sol, i, j)
                if delta < 0
                    sol.tour[i:j] = reverse(sol.tour[i:j])
                    sol.obj += delta
                    improved = true
                end
            end
        end
    end
    return nothing
end
