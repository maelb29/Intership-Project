# Dossier contenant les fonctions sur les variables globales pour le VRP
#
#
#include("Instances/inst_N4_W0106_M8_P40_Id2_TD30_DistOrig.jl") # pour les tests


# Fonction qui permet de créer la variable G associé au chemin S
function ForwardTimeSlackVRP(S, travelTime, TWB, TWE, meanTaskduration, shift, shiftB, shiftE, maxO)
    n = length(S)
    G1 = Float64[]  # Arrival times
    G2 = Float64[]  # Latest departure times
    G3 = Float64[]  # Waiting times

    # === G[1] : arrival times ===
    push!(G1, shiftB[shift])  # départ du dépôt

    for i in 2:n
        prev = S[i-1]
        curr = S[i]

        prev_duration = prev == 0 ? 0.0 : meanTaskduration[prev]
        window_start  = curr == 0 ? shiftB[shift] : TWB[curr]

        arrival = max(window_start, G1[i-1] + travelTime[prev+1, curr+1] + prev_duration)
        push!(G1, arrival)
    end

    # === G[2] : latest feasible times (retour arrière) ===
    resize!(G2, n)
    G2[n] = S[n] == 0 ? shiftE[shift] + maxO : TWE[S[n]]

    for i in (n-1):-1:1
        curr = S[i]
        nxt = S[i+1]

        curr_tw_end  = curr == 0 ? shiftE[shift] + maxO : TWE[curr]
        nxt_duration = nxt == 0 ? 0.0 : meanTaskduration[nxt]

        curr_duration = curr == 0 ? 0.0 : meanTaskduration[curr]
    G2[i] = min(curr_tw_end + curr_duration, G2[i+1] - travelTime[curr+1, nxt+1] - nxt_duration)
    end

    # === G[3] : waiting times ===
    push!(G3, 0.0)  # pas d'attente au départ

    for k in 2:n-1
        prev = S[k-1]
        curr = S[k]

        prev_duration = prev == 0 ? 0.0 : meanTaskduration[prev]
        window_start  = curr == 0 ? shiftB[shift] : TWB[curr]

        arrival    = G1[k-1] + prev_duration + travelTime[prev+1, curr+1]
        wait_time  = max(window_start - arrival, 0.0)
        push!(G3, wait_time)
    end

    push!(G3, 0.0)  # pas d’attente à l’arrivée dépôt retour

    return [G1, G2, G3]
end


# Fonction auxiliaire qui permet de trouver la clé d'un dictionnaire à partir de sa valeur (utilisé dans les fonction update_remove et update_add)
function trouver_patient(task, taskByPatient)
    for (patient, tasks) in taskByPatient
        if task in tasks
            return patient
        end
    end
    return nothing
end

# Fonction qui renvoie la nouvelle valeur de S et G après suppression d'un sommet et qui MAJ les variables telles que nurse_workload et nurse_patient
function update_remove(G, S, i, travelTime, TWB, TWE, meanTaskduration, shift, shiftB, shiftE, nurse, current_nurse_workload, nurse_by_patient, taskbypatient)

    v = S[i]
    # nouvelle duration potentiellement perdue
    duration = shiftE[shift] - shiftB[shift]
    # Suppression
    deleteat!(S, i)
    deleteat!(G[1], i)
    deleteat!(G[2], i)

    # === Recalcul G[1] : arrival times ===
    for k in i:length(S)
        prev = S[k-1]
        curr = S[k]
        prev_duration = prev == 0 ? 0.0 : meanTaskduration[prev]
        window_start = curr == 0 ? shiftB[shift] : TWB[curr]
        new_arrival = max(window_start, G[1][k-1] + travelTime[prev+1, curr+1] + prev_duration)

        if G[1][k] == new_arrival
            break
        end
        G[1][k] = new_arrival
    end

    # === Recalcul G[2] : latest times ===
    for k in (i-1):-1:1
        curr = S[k]
        nxt = S[k+1]
        curr_tw_end = curr == 0 ? shiftE[shift] : TWE[curr]
        nxt_duration = nxt == 0 ? 0.0 : meanTaskduration[nxt]
        new_latest = min(curr_tw_end, G[2][k+1] - travelTime[curr+1, nxt+1] - nxt_duration)

        if G[2][k] == new_latest
            break
        end
        G[2][k] = new_latest
    end

    # === Recalcul G[3] : waiting times ===
    G[3] = [0.0]
    for k in 2:length(S)-1
        prev = S[k-1]
        curr = S[k]
        prev_duration = prev == 0 ? 0.0 : meanTaskduration[prev]
        window_start = curr == 0 ? shiftB[shift] : TWB[curr]
        arrival = G[1][k-1] + prev_duration + travelTime[prev+1, curr+1]
        wait_time = max(window_start - arrival, 0.0)
        push!(G[3], wait_time)
    end
    push!(G[3], 0.0)

    # Cas où l'infirmière ne travaille plus sur cette journée donc perd sa shift 
    if length(S) == 2
        current_nurse_workload[nurse] = get(current_nurse_workload, nurse, 0.0) - duration
    end

    # update de nurse_by_patient
    patient = trouver_patient(v, taskbypatient)
    #deleatead!(nurse_by_patient[patient], findfirst(==(nurse), nurse_by_patient[patient]))

    # Vérifier si ce patient est bien suivi dans le dico
    if haskey(nurse_by_patient, patient) && haskey(nurse_by_patient[patient], nurse)
        if nurse_by_patient[patient][nurse] > 1
            nurse_by_patient[patient][nurse] -= 1
        else
            delete!(nurse_by_patient[patient], nurse)
        end
    end

    return S, G
end


# Fonction qui renvoie la nouvelle valeur de S et G après l'ajout d'un sommet et qui MAJ les variables telles que nurse_workload et nurse_patient
function update_add(G, S, v, i, travelTime, TWB, TWE, meanTaskduration, shiftB, shiftE, shift, nurse, nurse_by_patient, taskBypatient)

    # Nouveau temps de travail potentiel si la nurse n'avait pas de route
    duration = shiftE[shift] - shiftB[shift]

    # Insertion
    insert!(S, i, v)
    insert!(G[1], i, 0.0)
    insert!(G[2], i, 0.0)

    # === Update G[1] : heures d'arrivée ===
    for k in i:length(S)
        prev = S[k-1]
        curr = S[k]

        prev_duration = prev == 0 ? 0.0 : meanTaskduration[prev]
        window_start = curr == 0 ? shiftB[shift] : TWB[curr]

        arrival = max(window_start, G[1][k-1] + travelTime[prev+1, curr+1] + prev_duration)

        if G[1][k] == arrival
            break
        end
        G[1][k] = arrival
    end

    # === Update G[2] : dernières heures possibles ===
    curr = S[i]
    nxt = S[i+1]

    curr_tw_end = curr == 0 ? shiftE[shift] : TWE[curr]
    nxt_duration = nxt == 0 ? 0.0 : meanTaskduration[nxt]

    G[2][i] = min(curr_tw_end + (curr == 0 ? 0.0 : meanTaskduration[curr]),
                  G[2][i+1] - travelTime[curr+1, nxt+1] - nxt_duration)

    for k in (i-1):-1:1
        curr = S[k]
        nxt = S[k+1]

        curr_tw_end = curr == 0 ? shiftE[shift] : TWE[curr]
        nxt_duration = nxt == 0 ? 0.0 : meanTaskduration[nxt]

        new_depart = min(curr_tw_end + (curr == 0 ? 0.0 : meanTaskduration[curr]),
                         G[2][k+1] - travelTime[curr+1, nxt+1] - nxt_duration)

        if G[2][k] == new_depart
            break
        end
        G[2][k] = new_depart
    end

    # === Update G[3] : temps d'attente ===
    G[3] = [0.0]
    for k in 2:length(S)-1
        prev = S[k-1]
        curr = S[k]

        prev_duration = prev == 0 ? 0.0 : meanTaskduration[prev]
        window_start = curr == 0 ? shiftB[shift] : TWB[curr]

        arrival = G[1][k-1] + prev_duration + travelTime[prev+1, curr+1]
        wait_time = max(window_start - arrival, 0.0)
        push!(G[3], wait_time)
    end
    push!(G[3], 0.0)

    # update de nurse_by_patient
    patient = trouver_patient(v, taskBypatient)
    
    if !haskey(nurse_by_patient, patient)
        nurse_by_patient[patient] = Dict{Int,Int}()
    end

    if haskey(nurse_by_patient[patient], nurse)
        nurse_by_patient[patient][nurse] += 1
    else
        nurse_by_patient[patient][nurse] = 1
    end

    # Ici on a pas besoin de mettre à jour les heures de travail des infirmière car les infirmières ne sont pas changées de shift
    
    return S, G
end


# Fonction qui vérifie si l'insertion d'une tâche v à ième position du chemin S est possible
# Si l'insertion est possible renvoie true et le coût de cette insertion (utile pour les opérateurs de construction par la suite) sinon renvoie false
function FeasibleInsertion(S, G, travelTime, v, i, TWB, TWE, shiftB, meanTaskduration, shift, external_cost, travel_cost, day, nurse)

    cur = S[i]
    next = S[i+1]

    if 2 ≤ i ≤ length(S)-2
        w = max(TWB[v], G[1][i] + travelTime[cur+1,v+1] + meanTaskduration[cur])
        delta = max(TWB[next], w + meanTaskduration[v] + travelTime[v+1,next+1])

        if (w ≤ TWE[v]) && (delta + meanTaskduration[next] ≤ G[2][i+1])
            cost_increase = (travelTime[cur+1,v+1] + travelTime[v+1,next+1] - travelTime[cur+1,next+1])*travel_cost[nurse, day] - external_cost[v] 
            return true, cost_increase
        end

    elseif i == 1
        w = max(TWB[v], G[1][i] + travelTime[cur+1,v+1])
        delta = max(TWB[next], w + meanTaskduration[v] + travelTime[v+1,next+1])

        if (w ≤ TWE[v]) && (delta + meanTaskduration[next] ≤ G[2][i+1])
            cost_increase = (travelTime[cur+1,v+1] + travelTime[v+1,next+1] - travelTime[cur+1,next+1])*travel_cost[nurse, day] - external_cost[v] 
            return true, cost_increase
        end

    else # i == length(S)-1
        w = max(TWB[v], G[1][i] + travelTime[cur+1,v+1] + meanTaskduration[cur])
        delta = max(shiftB[shift], w + meanTaskduration[v] + travelTime[v+1,next+1])

        if (w ≤ TWE[v]) && (delta ≤ G[2][i+1])
        
            cost_increase = (travelTime[cur+1,v+1] + travelTime[v+1,next+1] - travelTime[cur+1,next+1])*travel_cost[nurse, day] - external_cost[v] 
            return true, cost_increase
        end
    end

    return false, Inf
end


