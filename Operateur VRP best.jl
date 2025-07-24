using Random
include("Variable Globales VRP best.jl") 
# include("Instances_generees/inst_N4_W0106_M8_P40_Id2_TD30_DistOrig.jl") # sert uniquement pour les tests


## OPERATEUR DE DESTRUCTION


function randomRemove(Planning, request_banque, q, travelTime, TWB, TWE, meanTaskduration, 
    shiftB, shiftE, current_nurse_workload, nurse_by_patient, taskbypatient)

    removed_total = 0
    max_attempts = 500
    attempt = 0

    while removed_total < q && attempt < max_attempts
        attempt += 1
        # Choisir un jour aléatoire
        jour_idx = rand(1:length(Planning))
        PlanningJour = Planning[jour_idx]
        request_banque_jour = request_banque[jour_idx]

        # Filtrer les tournées valides
        valid_routes = [i for i in 1:length(PlanningJour) if !ismissing(PlanningJour[i][3]) && length(PlanningJour[i][3]) > 2]

        if isempty(valid_routes)
            continue
        end

        # Sélectionner une tournée aléatoire
        idx = rand(valid_routes)
        path = PlanningJour[idx][3]
        G = PlanningJour[idx][4]
        shift = PlanningJour[idx][2]
        nurse = PlanningJour[idx][1]

        # Choisir une tâche à retirer (hors dépôts)
        removable_indices = 2:(length(path) - 1)
        task_index = rand(removable_indices)
        removed_task = path[task_index]

        # Mise à jour du chemin et du G, avec MAJ des heures
        new_path, new_G = update_remove(
            G, path, task_index,
            travelTime, TWB, TWE, meanTaskduration,
            shift, shiftB, shiftE,
            nurse, current_nurse_workload, nurse_by_patient, taskbypatient
        )

        # Mettre à jour la tournée dans le planning
        PlanningJour[idx][3] = new_path
        PlanningJour[idx][4] = new_G

        if length(new_path) == 2
            PlanningJour[idx][2] = missing
        end

        # Stocker la tâche retirée
        push!(request_banque_jour, removed_task)

        removed_total += 1

        # Réinjecter les structures modifiées
        Planning[jour_idx] = PlanningJour
        request_banque[jour_idx] = request_banque_jour
    end
end



# Fonction qui calcule le cout de retrait d'un sommet à la position i dans un chemin
function computeRemoveIncrease(S, i, travelTime, external_cost, nurse, hire_cost, travel_cost, day, shift)

    prev, curr, nxt = S[i-1], S[i], S[i+1]
    cost_increase = external_cost[curr] + travelTime[prev+1, nxt+1]*travel_cost[nurse, day] -
                    (travelTime[prev+1, curr+1]*travel_cost[nurse, day] + travelTime[curr+1, nxt+1]*travel_cost[nurse, day]) # en approximation on ne prend pas en compte le temps d'overtime car ça demande de recalculer toute la solution

    #println("La shift du compute remove Increase : ", shift)
    if length(S) == 3
        cost_increase -= hire_cost[nurse, shift+1, day]
    end


    return cost_increase
end



# Opérateur worstRemoval
function worstRemoval(Planning, request_banque, q, travelTime, TWB, TWE,
     meanTaskduration, shiftB, shiftE, external_cost, hire_cost, travel_cost, p, taskByPatient, current_nurse_workload, nurse_by_patient)

    removed_total = 0

    while removed_total < q
        # === Construction des candidats dans tout le planning ===
        candidates = []

        for day in 1:length(Planning)
            PlanningJour = Planning[day]

            for route_idx in 1:length(PlanningJour)
                nurse = PlanningJour[route_idx][1]
                shift = PlanningJour[route_idx][2]
                path = PlanningJour[route_idx][3]

                if shift === missing || length(path) ≤ 2
                    continue
                end
            
                #println()
                #println("la shift avant de rentrer dans compute remove increase : ", shift)
                for i in 2:length(path)-1
                    cost = computeRemoveIncrease(path, i, travelTime, external_cost, nurse, hire_cost, travel_cost, day, shift)
                    push!(candidates, (cost, day, route_idx, i))
                end
            end
        end

        # Si aucun candidat disponible, on arrête
        if isempty(candidates)
            break
        end

        # === Tri des candidats et sélection biaisée
        sort!(candidates, by = x -> x[1])
        y = rand()
        idx = Int(floor(y^p * length(candidates))) + 1
        cost, day, route_idx, task_idx = candidates[idx]

        # === Suppression du sommet sélectionné
        PlanningJour = Planning[day]
        path = PlanningJour[route_idx][3]
        G = PlanningJour[route_idx][4]
        shift = PlanningJour[route_idx][2]
        task_removed = path[task_idx]
        nurse = PlanningJour[route_idx][1]

        new_path, new_G = update_remove(G, path, task_idx, travelTime, TWB, TWE, meanTaskduration, shift, shiftB,
         shiftE, nurse, current_nurse_workload, nurse_by_patient, taskByPatient)

        PlanningJour[route_idx][3] = new_path
        PlanningJour[route_idx][4] = new_G

        if length(new_path) ≤ 2
            PlanningJour[route_idx][2] = missing
        end

        # === Mise à jour de la banque de requêtes
        push!(request_banque[day], task_removed)
        removed_total += 1

        # Réinjecter le jour mis à jour
        Planning[day] = PlanningJour
    end
end



# Opérateur waitingTimeRemoval
function waitingTimeRemoval(
    Planning, request_banque, travelTime, TWB, TWE, meanTaskduration,
    shiftB, shiftE, q, pworst, nurse_workload, nurse_by_patient, taskbypatient)

    removed = Int[]

    while length(removed) < q
        candidates = []

        for day in 1:length(Planning)
            PlanningJour = Planning[day]

            for route_idx in 1:length(PlanningJour)
                shift = PlanningJour[route_idx][2]
                path = PlanningJour[route_idx][3]
                G = PlanningJour[route_idx][4]

                if shift === missing || length(path) ≤ 2
                    continue
                end

                waiting_times = G[3]
                for i in 2:length(path)-1
                    wait = waiting_times[i]
                    push!(candidates, (wait, day, route_idx, i))
                end
            end
        end

        if isempty(candidates)
            break
        end

        # Trier par temps d’attente décroissant
        sort!(candidates, by = x -> -x[1])

        # Sélection biaisée (pworst)
        y = rand()
        idx = clamp(Int(floor(y^pworst * length(candidates))) + 1, 1, length(candidates))
        _, day, route_idx, task_idx = candidates[idx]

        # === Suppression du sommet
        PlanningJour = Planning[day]
        path = PlanningJour[route_idx][3]
        G = PlanningJour[route_idx][4]
        shift = PlanningJour[route_idx][2]
        nurse = PlanningJour[route_idx][1]

        task_removed = path[task_idx]

        new_path, new_G = update_remove(G, path, task_idx, travelTime, TWB, TWE, meanTaskduration, shift, shiftB, 
        shiftE, nurse, nurse_workload, nurse_by_patient, taskbypatient)

        Planning[day][route_idx][3] = new_path
        Planning[day][route_idx][4] = new_G

        if length(new_path) ≤ 2
            Planning[day][route_idx][2] = missing
        end

        push!(request_banque[day], task_removed)
        push!(removed, task_removed)
    end
end




# Fonction auxiliaire pour calculer la distance de relatedness entre deux tâches pour le shawRemoval
function relatedness(i, j, TWB, travelTime, Skills)
    # Distance de routing (travel time)
    dist = travelTime[i+1, j+1]

    # Différence de temps de début de fenêtre temporelle
    time_diff = abs(TWB[i] - TWB[j])

    # Distance de compétence entre tâches :
    # On définit la distance de Hamming comme suit :
    # Si deux tâches partagent des infirmières en commun, elles sont plus proches.
    nurses_i = [nurse for (nurse, tasks) in Skills if i in tasks]
    nurses_j = [nurse for (nurse, tasks) in Skills if j in tasks]
    #println("Ensemble nurses_i : ", nurses_i)
    #println("Ensemble nurses_j : ", nurses_j)

    # Calcul de la symdiff manuelle
    skill_distance = 0
    for n in nurses_i
        if !(n in nurses_j)
            skill_distance += 1
        end
    end
    for n in nurses_j
        if !(n in nurses_i)
            skill_distance += 1
        end
    end
    #println("Symdiff : ", skill_distance)

    # Poids (modifiables)
    δ = 1.0     # poids de la distance
    φ = 1.0     # poids différence TW
    μ = 1.0     # poids compétence

    return δ * dist + φ * time_diff + μ * skill_distance
end



# Opérateur relatedRemoval (il s'agit du ShawRemoval)
function relatedRemoval(Planning, request_banque, q, TWB, TWE, travelTime, Skills,
    prelated, meanTaskduration, shiftB, shiftE, nurse_workload, nurse_by_patient, taskbypatient)

    removed_ids = Int[]
    to_remove = []  # Stocke les quadruplets (day, route_idx, task_idx, task_id)
    remaining_tasks = []

    # === Étape 0 : collecter toutes les tâches de tous les jours
    for day in 1:length(Planning)
        for (route_idx, nurse) in enumerate(Planning[day])
            path = nurse[3]
            if !ismissing(path) && length(path) > 2
                for i in 2:length(path)-1
                    task = path[i]
                    push!(remaining_tasks, (day, route_idx, i, task))
                end
            end
        end
    end

    if isempty(remaining_tasks)
        return removed_ids
    end

    # === Étape 1 : choisir une tâche initiale aléatoirement
    first = rand(remaining_tasks)
    push!(to_remove, first)
    push!(removed_ids, first[4])  # task_id
    deleteat!(remaining_tasks, findfirst(==(first), remaining_tasks))

    # === Étape 2 : suppressions guidées par similarité
    while length(removed_ids) < q && !isempty(remaining_tasks)
        anchor = rand(removed_ids)

        rel_list = [(relatedness(anchor, t[4], TWB, travelTime, Skills), t) for t in remaining_tasks]
        sort!(rel_list, by = x -> x[1])

        y = rand()
        idx = clamp(Int(floor(y^prelated * length(rel_list))) + 1, 1, length(rel_list))
        _, candidate = rel_list[idx]

        push!(to_remove, candidate)
        push!(removed_ids, candidate[4])
        deleteat!(remaining_tasks, findfirst(==(candidate), remaining_tasks))
    end

    # === Étape 3 : appliquer les suppressions triées par indices décroissants
    grouped = Dict{Tuple{Int, Int}, Vector{Tuple{Int, Int}}}()  # (day, route_idx) => [(task_idx, task_id)]

    for (day, route_idx, task_idx, task_id) in to_remove
        key = (day, route_idx)
        grouped[key] = get(grouped, key, Tuple{Int, Int}[])
        push!(grouped[key], (task_idx, task_id))
    end

    for ((day, route_idx), entries) in grouped
        sort!(entries, by = x -> -x[1])  # supprimer du plus grand index vers le plus petit
        for (task_idx, task_id) in entries
            path  = Planning[day][route_idx][3]
            G     = Planning[day][route_idx][4]
            shift = Planning[day][route_idx][2]
            nurse = Planning[day][route_idx][1]

            new_path, new_G = update_remove(G, path, task_idx, travelTime, TWB, TWE, meanTaskduration, shift, shiftB, shiftE, 
            nurse, nurse_workload, nurse_by_patient, taskbypatient)

            Planning[day][route_idx][3] = new_path
            Planning[day][route_idx][4] = new_G
            if length(new_path) <= 2
                Planning[day][route_idx][2] = missing
            end
            push!(request_banque[day], task_id)
        end
    end
end



# Fonction qui renvoie un tableau qui décrit le nombre d'infirmière qui s'occupe d'un patient (utilisé pour l'opérateur concentratedRemoval)
function tableau_affectation(Planning, taskByPatient, P)
    D = length(Planning)
    tableau = Array{Int64}(undef, P, D + 1)
    fill!(tableau, 0)

    # Pour chaque jour
    for d in 1:D
        jour = Planning[d]
        # Pour chaque tournée d’infirmière ce jour-là
        for route in jour
            nurse = route[1]
            tasks = route[3]

            if ismissing(tasks)
                continue
            else
                for task in tasks
                    for (patient, patient_tasks) in taskByPatient
                        #println(task, " ", patient, " ", patient_tasks)
                        if task in patient_tasks
                            tableau[patient, d] = nurse
                        end
                    end
                end
            end
        end
    end

    # Dernière colonne : nombre total d'infirmières différentes par patient
    for patient in 1:P
        nurses = Set{Int}()
        for d in 1:D
            nurse = tableau[patient, d]
            if nurse != 0
                push!(nurses, nurse)
            end
        end
        tableau[patient, D+1] = length(nurses)
    end

    return tableau
end



# Opérateur de destruction qui supprime les taches des patients qui sont assignés au plus d'infirmières
# Attention q ici ne correspond pas au nombre de tache supprimée mais au nombre de patient 
function concentratedRemoval(
    Planning, request_banque, q, travelTime, TWB, TWE, meanTaskduration,
    shiftB, shiftE, taskByPatient, current_nurse_workload, P, nurse_by_patient)

    # 1. Construire le tableau de présence
    affectation = tableau_affectation(Planning, taskByPatient, P)

    # 2. Identifier les patients les plus "dispersés"
    scores = [(p, affectation[p, end]) for p in 1:P]  # patient => nb d'infirmières
    sort!(scores, by = x -> -x[2])  # décroissant

    patients_cibles = [scores[i][1] for i in 1:min(q, length(scores))]

    # 3. Collecte des tâches à retirer pour ces patients
    tasks_to_remove = Set{Int}()
    for p in patients_cibles
        for task in taskByPatient[p]
            push!(tasks_to_remove, task)
        end
    end

    # 4. Suppression des tâches dans le planning
    for day in 1:length(Planning)
        for route_idx in 1:length(Planning[day])
            nurse = Planning[day][route_idx][1]
            shift = Planning[day][route_idx][2]
            path = Planning[day][route_idx][3]
            G = Planning[day][route_idx][4]


            # Vérifie que path n'est pas missing
            if ismissing(path)
                continue
            end


            to_remove = [i for i in 2:length(path)-1 if path[i] in tasks_to_remove]

            # Supprimer en partant du fond
            for idx in reverse(to_remove)
                task_removed = path[idx]
                path, G = update_remove(G, path, idx, travelTime, TWB, TWE,
                meanTaskduration, shift, shiftB, shiftE, nurse, current_nurse_workload, nurse_by_patient, taskByPatient)

                push!(request_banque[day], task_removed)
            end

            # Mettre à jour la route
            Planning[day][route_idx][3] = path
            Planning[day][route_idx][4] = G
            if length(path) ≤ 2
                Planning[day][route_idx][2] = missing
            end
        end
    end
end





## OPERATEUR DE CONSTRUCTION

# Fonction qui calcule le cout d'insertion d'une tache à la position i dans le chemin S pour length(S)>2
function ComputeCostIncrease(S, task, i, travelTime, external_cost, nurse, day, travel_cost)
    total_cost = 0
    # Pas besoin de prendre le hiring cost ici puisque le chemin est déjà existant donc il n'y pas de cout d'embauche car l'infirmière travaille déjà
    total_cost += (travelTime[S[i-1],task] + travelTime[task,S[i]] - travelTime[S[i-1],S[i]])*travel_cost[nurse, day] - external_cost[task]  

    return total_cost
end



# Fonction qui effectue l'insertion d'un sommet v dans un chemin qui a pour shift missing
# Elle remplace updateG_add 
function EmptyInsertion(task, travelTime, TWB, TWE, meanTaskduration, shift, shiftB, shiftE, maxO)
    S_new = [0,task,0]
    G_new = ForwardTimeSlackVRP(S_new, travelTime, TWB, TWE, meanTaskduration, shift, shiftB, shiftE, maxO) 

    #println("valeur de la shift : ", shift)
    #println("La valeur de S_new : ", S_new)
    #println("La valeur de G : ", G)
    return S_new, G_new
end


# Fonction qui regarde si la création d'une shift pour la nurse au jour day est compatible pour les temps de repos de la nurse
function RestCheckEmpty(Planning, shift, nurse, shiftE, shiftB, day, rest_hour)
    if day == 1
        shift_next = Planning[day+1][nurse][2]

        if ismissing(shift_next)
            return true
        else
            gap_shift_next = shiftB[shift_next] - shiftE[shift] + 24*60
            return gap_shift_next >= rest_hour * 60
        end

    elseif day == length(Planning)
        shift_prev = Planning[day-1][nurse][2]

        if ismissing(shift_prev)
            return true
        else
            gap_shift_prev = shiftB[shift] - shiftE[shift_prev] + 24*60
            return gap_shift_prev >= rest_hour * 60
        end

    else
        shift_next = Planning[day+1][nurse][2]
        shift_prev = Planning[day-1][nurse][2]

        ok_next = true
        ok_prev = true

        if !ismissing(shift_next)
            gap_shift_next = shiftB[shift_next] - shiftE[shift] + 24*60
            ok_next = gap_shift_next >= rest_hour * 60
        end

        if !ismissing(shift_prev)
            gap_shift_prev = shiftB[shift] - shiftE[shift_prev] + 24*60
            ok_prev = gap_shift_prev >= rest_hour * 60
        end

        return ok_next && ok_prev
    end
end



# Fonction qui vérifie si l'insertion de la tache est possible par rapport à aux fenêtres horaires de la tâche et de la shift
function windowInsertion(task, travelTime, TWB, TWE, shiftB, shiftE, shift, meanTaskduration)
    # Vérification initiale : chevauchement fenêtre temporelle / shift
    if TWB[task] <= shiftE[shift] && TWE[task] >= shiftB[shift]
        # Temps d'arrivée à la tâche en partant du dépôt
        arrival_task = max(shiftB[shift] + travelTime[1, task+1], TWB[task])
        #println("arrival time vaut : ", arrival_task)
        if arrival_task <= TWE[task]
            # Heure retour au dépôt après la tâche
            arrival_depot = arrival_task + meanTaskduration[task] + travelTime[task+1,1]

            # On vérifie que tout rentre dans la shift
            return arrival_depot <= shiftE[shift]
        else
            return false
        end
    else
        return false
    end
end


#=
# Fonction qui calcule la meilleur shift en calculant le nombre de tache potentiel que la shift peut recueillir en plus
function IdealShiftAssignment(Planning, task, banque_jour, TWB, TWE, shiftB, shiftE, skills, meanTaskduration, nurse, day, rest_hour, max_weekly_work, nurse_workload)
    best_shift = -1
    best_score = -1

    for k in 1:length(shiftB)
        # Est-ce que la shift est compatible avec la tâche principale ?
        if windowInsertion(task, travelTime, TWB, TWE, shiftB, shiftE, k, meanTaskduration) &&
            RestCheckEmpty(Planning, k, nurse, shiftE, shiftB, day, rest_hour)

            augmentation_duration = shiftE[k] - shiftB[k]
            #println("Ca passe pour la shift ", k)

            if nurse_workload[nurse] + augmentation_duration <= max_weekly_work
                score = 0
                #println("Ca passe encore pour la shift ", k)

                # On regarde combien d'autres tâches de la request banque cette shift peut potentiellement accueillir
                for other_task in banque_jour
                    if TWB[other_task] <= shiftE[k] && TWE[other_task] >= shiftB[k] && (other_task in skills[nurse])
                        score += 1
                        #println("On a ajouté 1 au score pour la shift ", k, " et la tâche ", other_task)
                    end
                end

                #println("score pour la shift ", k, "  ", score)
                #println()

                # Mise à jour du meilleur choix
                if score > best_score
                    best_score = score
                    best_shift = k
                end
            end
        end
    end

    return best_shift
end
=#



# Fonction qui calcule la meilleur shift en calculant le nombre de tache potentiel que la shift peut recueillir en plus
function IdealShiftAssignment(Planning, task, TWB, TWE, shiftB, shiftE, meanTaskduration, nurse, day, rest_hour, max_weekly_work, nurse_workload)
    possible_shift = []

    for k in 1:length(shiftB)
        # Est-ce que la shift est compatible avec la tâche principale ?
        if windowInsertion(task, travelTime, TWB, TWE, shiftB, shiftE, k, meanTaskduration) &&
            RestCheckEmpty(Planning, k, nurse, shiftE, shiftB, day, rest_hour)

            augmentation_duration = shiftE[k] - shiftB[k]
            #println("Ca passe pour la shift ", k)

            if nurse_workload[nurse] + augmentation_duration <= max_weekly_work
                push!(possible_shift, k)
            end

        end
    end

    #println("L'ensemble des shifts possibles est : ", possible_shift)

    if length(possible_shift) > 0
        random_shift = rand(possible_shift)
        return random_shift
    else
        return -1
    end
end




# Opérateur greedyInsertion
function GreedyInsertion(
    Planning, request_banque, travelTime, TWB, TWE,
    meanTaskduration, shiftB, shiftE, external_cost, hire_cost, travel_cost,
    taskByPatient, maxO, skills, nurse_workload, max_weekly_work, rest_hour, nurse_by_patient, maxNumberNurse)

    while any(!isempty, request_banque)
        best_cost = Inf
        best_task = -1
        best_day = -1
        best_route_idx = -1
        best_insert_idx = -1
        best_shift = -1

        for d in 1:length(Planning)
            PlanningJour = Planning[d]
            banque_jour = request_banque[d]

            for task in banque_jour
                for route_idx in 1:length(PlanningJour)
                    nurse = PlanningJour[route_idx][1]
                    shift = PlanningJour[route_idx][2]
                    path = PlanningJour[route_idx][3]
                    G = PlanningJour[route_idx][4]

                    # On vérifie l'infirmière possède bien la skills
                    if !(task in skills[nurse])
                        continue
                    end

                    # On vérifie que le nombre d'infirmière associé au patient ne dépasse pas une certaine valeur
                    patient = trouver_patient(task, taskByPatient)
                    if haskey(nurse_by_patient, patient)
                        # Si le nurse n'est pas déjà affecté à ce patient, on vérifie la limite
                        if !(haskey(nurse_by_patient[patient], nurse)) && length(nurse_by_patient[patient]) >= maxNumberNurse[patient] 
                            continue  # on saute cette affectation
                        end
                    end

                   

                    if shift === missing
                        current_workload = get(nurse_workload, nurse, 0.0)
                        # Cas tournée vide → on tente une insertion directe
                        s = IdealShiftAssignment(Planning, task, TWB, TWE, shiftB, shiftE, meanTaskduration,
                                                nurse, d, rest_hour, max_weekly_work, nurse_workload)
                        if s !== -1
                            cost = (travelTime[1, task+1] + travelTime[task+1, 1])*travel_cost[nurse, d] - external_cost[task] +
                                     + hire_cost[nurse, s+1, d]

                            if cost < best_cost
                                best_cost = cost
                                best_task = task
                                best_day = d
                                best_route_idx = route_idx
                                best_insert_idx = :empty
                                best_shift = s
                            end
                        end

                    else
                        for i in 1:length(path)-1
                            feasible, cost = FeasibleInsertion(
                                path, G, travelTime, task, i,
                                TWB, TWE, shiftB, meanTaskduration, shift,
                                external_cost, travel_cost, d, nurse,
                            )

                            if feasible && cost < best_cost
                                best_cost = cost
                                best_task = task
                                best_route_idx = route_idx
                                best_insert_idx = i
                                best_shift = shift
                                best_day = d
                            end
                        end
                    end
                end
            end
        end

        if best_task === -1
            break  # Rien d'insérable
        end

        # Mise à jour du planning
        route = Planning[best_day][best_route_idx]
        nurse = route[1]

        if best_insert_idx == :empty
            S_new, G_new = EmptyInsertion(best_task, travelTime, TWB, TWE, meanTaskduration, best_shift, shiftB, shiftE, maxO)
            route[3] = S_new
            route[4] = G_new
            route[2] = best_shift

            Planning[best_day][best_route_idx] = route

            # Mise à jour de la durée de la nouvelle tournée
            duration = shiftE[best_shift] - shiftB[best_shift]
            nurse_workload[nurse] = get(nurse_workload, nurse, 0.0) + duration

            # update de nurse_by_patient
            patient = trouver_patient(best_task, taskByPatient)

            if !haskey(nurse_by_patient, patient)
                nurse_by_patient[patient] = Dict{Int,Int}()
            end

            if haskey(nurse_by_patient[patient], nurse)
                nurse_by_patient[patient][nurse] += 1
            else
                nurse_by_patient[patient][nurse] = 1
            end

            #println("La tache ", best_task, " a été ajouté au chemin vide ", best_route_idx, " le jour ", best_day)
        else
            S_new, G_new = update_add(
                route[4], route[3], best_task, best_insert_idx + 1,
                travelTime, TWB, TWE, meanTaskduration, shiftB, shiftE, best_shift, nurse, nurse_by_patient, taskByPatient
            )

            Planning[best_day][best_route_idx][3] = S_new
            Planning[best_day][best_route_idx][4] = G_new

            #println("La tache ", best_task, " a été ajouté au chemin ", best_route_idx, " le jour ", best_day)

        end


        deleteat!(request_banque[best_day], findfirst(==(best_task), request_banque[best_day]))

        #println("Planning après ajout :", Planning)
        #println()
    end
end



# Opérateur k-regret
function RegretInsertion(
    Planning, request_banque, travelTime, TWB, TWE,
    meanTaskduration, shiftB, shiftE, external_cost, hire_cost, travel_cost,
    taskByPatient, maxO, skills, nurse_workload, max_weekly_work, k, rest_hour, nurse_by_patient, maxNumberNurse)
    

    while any(!isempty, request_banque)
        best_regret = -Inf
        best_task = -1
        best_day = -1
        best_choice = -1  # (cost, route_idx, insert_idx, shift)

        for d in 1:length(Planning)
            PlanningJour = Planning[d]
            banque_jour = request_banque[d]

            for task in banque_jour
                insertions = []

                for route_idx in 1:length(PlanningJour)
                    nurse = PlanningJour[route_idx][1]
                    shift = PlanningJour[route_idx][2]
                    path = PlanningJour[route_idx][3]
                    G = PlanningJour[route_idx][4]

                    if !(task in skills[nurse])
                        continue
                    end

                    # On vérifie que le nombre d'infirmière associé au patient ne dépasse pas une certaine valeur
                    patient = trouver_patient(task, taskByPatient)
                    if haskey(nurse_by_patient, patient)
                        # Si le nurse n'est pas déjà affecté à ce patient, on vérifie la limite
                        if !(haskey(nurse_by_patient[patient], nurse)) && length(nurse_by_patient[patient]) >= maxNumberNurse[patient]
                            continue  # on saute cette affectation
                        end
                    end



                    if shift === missing
                        s = IdealShiftAssignment(Planning, task, TWB, TWE, shiftB, shiftE, meanTaskduration, nurse, d,
                                                rest_hour, max_weekly_work, nurse_workload)
                        if s !== -1
                            cost = (travelTime[1, task+1] + travelTime[task+1, 1])*travel_cost[nurse, d] - external_cost[task] + hire_cost[nurse, s+1, d]
                            push!(insertions, (cost, route_idx, :empty, s))
                        end

                    else
                        for i in 1:length(path)-1
                            feasible, cost = FeasibleInsertion(
                                path, G, travelTime, task, i,
                                TWB, TWE, shiftB, meanTaskduration, shift,
                                external_cost, travel_cost, d, nurse
                            )
                            if feasible
                                push!(insertions, (cost, route_idx, i, shift))
                            end
                        end
                    end
                end

                if isempty(insertions)
                    continue
                end

                sorted_insertions = sort(insertions, by = x -> x[1])

                regret = if length(sorted_insertions) == 1
                    Inf
                else
                    n = min(k, length(sorted_insertions))
                    sum(sorted_insertions[k][1] for k in 2:n) - (n - 1) * sorted_insertions[1][1]
                end

                if regret > best_regret
                    best_regret = regret
                    best_task = task
                    best_choice = sorted_insertions[1]
                    best_day = d
                end
            end
        end

        if best_task === -1
            break  # Rien d’insérable
        end

        cost, route_idx, insert_idx, shift = best_choice
        route = Planning[best_day][route_idx]
        nurse = route[1]

        if insert_idx == :empty
            S_new, G_new = EmptyInsertion(best_task, travelTime, TWB, TWE, meanTaskduration, shift, shiftB, shiftE, maxO)
            route[3] = S_new
            route[4] = G_new
            route[2] = shift
            Planning[best_day][route_idx] = route

            # update de nurse_workload
            augmentation_duration = shiftE[shift] - shiftB[shift]
            nurse_workload[nurse] = get(nurse_workload, nurse, 0.0) + augmentation_duration

            # update de nurse_by_patient
            patient = trouver_patient(best_task, taskByPatient)

            if !haskey(nurse_by_patient, patient)
                nurse_by_patient[patient] = Dict{Int,Int}()
            end

            if haskey(nurse_by_patient[patient], nurse)
                nurse_by_patient[patient][nurse] += 1
            else
                nurse_by_patient[patient][nurse] = 1
            end

        else
            S_new, G_new = update_add(
                route[4], route[3], best_task, insert_idx + 1,
                travelTime, TWB, TWE, meanTaskduration, shiftB, shiftE, shift, nurse, nurse_by_patient, taskByPatient
            )
            Planning[best_day][route_idx][3] = S_new
            Planning[best_day][route_idx][4] = G_new

        end

        deleteat!(request_banque[best_day], findfirst(==(best_task), request_banque[best_day]))
    end
end




