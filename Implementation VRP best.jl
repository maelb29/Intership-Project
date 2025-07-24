using Random 
using Plots
using LinearAlgebra
using Distances
using MultivariateStats
using Statistics
using BenchmarkTools


include("Variable Globales VRP best.jl")
include("Operateur VRP best.jl")
#include("Instances_generees/inst_N4_W0106_M5_P40_Id0_TD30_DistOrig.jl")
include("Instances_generees/inst_N4_W0106_M5_P40_Id2_TD30_DistOrig.jl")
#include("Instances_generees/inst_N5_W0106_M8_P60_Id2_TD30_DistOrig.jl")
#include("Instances_generees/inst_N11_W0106_M10_P128_Id0_TD30_DistOrig.jl")


# Fonction qui crée une solution initiale vide (Planning vide et request_banque pleine)
function InitialSolution(N, D, J, lastTaskonday)
    Planning = []
    request_banque = [] # Il s'agit ici de la request banque initiale

    # Remplissage de Planning
    for k in 1:D
        day = []
        for i in 1:N
            nurse = Any[i]
            push!(nurse, missing)
            push!(nurse, missing)
            push!(nurse, missing)
            push!(day, nurse)
        end
        push!(Planning, day)

        # Remplissage de la request_banque
        request_day = []

        if k == 1 
            for i in 1:lastTaskonday[k]
                push!(request_day, i)
            end
        elseif k < D  
            for i in lastTaskonday[k-1]+1:lastTaskonday[k]
                push!(request_day, i)
            end
        else
            for i in lastTaskonday[k-1]+1:J
                push!(request_day, i)
            end
        end

        push!(request_banque, request_day)

    end

    return Planning, request_banque
end


# Fonction qui calcule le cout par jour du VRP
function compute_VRP_solution(PlanningJour, request_banque_jour, travelTime, external_cost, hire_cost, overtime_cost, travel_cost, day, shiftB)
    cost = 0
    # Cout lié au chemin des infirmière
    for route in PlanningJour
        nurse = route[1]
        shift = route[2]
        if ismissing(shift)
            continue
        else
            path = route[3]
            G = [4]

            path_cost = 0


            # Cout lié à la durée du chemin
            for k in 1:length(path)-1
                path_cost += travelTime[path[k]+1,path[k+1]+1] * travel_cost[nurse, day]
            end

            # cout lié à l'assignation d'une shift à une infirmière
            path_cost += hire_cost[nurse, shift+1, day]

            # cout de pénalité si on fait de l'overtime 
            delta = G[1][end] - shiftB[shift]
            if delta > 0
                path_cost += delta*overtime_cost[nurse]
            end

            cost += path_cost
        end
    end

    # Cout lié à l'externalisation 
    for task in request_banque_jour
        external = external_cost[task]
        cost += external 
    end

    # cout lié à la journée de travail d'une infirmière


    return cost
end


# Fonction qui calcule le coût de la solution total
function compute_cost(Planning, request_banque, travelTime, external_cost, hire_cost, overtime_cost, travel_cost, shiftB)
    total_cost = 0
    for k in 1:length(Planning)
        cost = compute_VRP_solution(Planning[k], request_banque[k], travelTime, external_cost, hire_cost, overtime_cost, travel_cost, k, shiftB)

        total_cost += cost
    end

    return total_cost 
end


# Fonction qui calcule et affiche le détail des coûts pour chaque jour
function compute_VRP_solutionbis(PlanningJour, request_banque_jour, travelTime, external_cost, hire_cost, overtime_cost, travel_cost, day, taskByPatient, shiftB)
    total_cost = 0.0
    travel_cost_total = 0.0
    externalization_cost_total = 0.0
    hiring_cost_total = 0.0
    overtime_cost_total = 0.0

    for route in PlanningJour
        nurse = route[1]
        shift = route[2]
        if ismissing(shift)
            continue
        else
            path = route[3]
            G = [4]  # ← Je suppose que c'est temporaire, sinon il faut corriger ça

            path_cost = 0.0

            for k in 1:length(path)-1
                travel_c = travelTime[path[k]+1, path[k+1]+1] * travel_cost[nurse, day]
                path_cost += travel_c
                travel_cost_total += travel_c
            end

            hire_c = hire_cost[nurse, shift+1, day]
            path_cost += hire_c
            hiring_cost_total += hire_c

            delta = G[1][end] - shiftB[shift]
            if delta > 0
                overtime_c = delta * overtime_cost[nurse]
                path_cost += overtime_c
                overtime_cost_total += overtime_c
            end

            total_cost += path_cost
        end
    end

    # Externalisation
    for task in request_banque_jour
        external_c = external_cost[task]
        total_cost += external_c
        externalization_cost_total += external_c
    end

    return (
        total_cost = total_cost,
        travel_cost = travel_cost_total,
        externalization_cost = externalization_cost_total,
        hiring_cost = hiring_cost_total,
        overtime_cost = overtime_cost_total,
    )
end

# Foncion qui calcule le détaille complet des différents coûts pour 
function compute_costbis(Planning, request_banque, travelTime, external_cost, hire_cost, overtime_cost, travel_cost, taskByPatient, shiftB)
    total_cost = 0.0
    travel_cost_total = 0.0
    externalization_cost_total = 0.0
    hiring_cost_total = 0.0
    
    overtime_cost_total = 0.0

    for k in 1:length(Planning)
        result = compute_VRP_solutionbis(Planning[k], request_banque[k], travelTime, external_cost, hire_cost, overtime_cost, travel_cost, k, taskByPatient, shiftB)

        total_cost += result.total_cost
        travel_cost_total += result.travel_cost
        externalization_cost_total += result.externalization_cost
        hiring_cost_total += result.hiring_cost
        overtime_cost_total += result.overtime_cost
    end

    return (
        total_cost = total_cost,
        travel_cost = travel_cost_total,
        externalization_cost = externalization_cost_total,
        hiring_cost = hiring_cost_total,
        overtime_cost = overtime_cost_total,
    )
end


# Fonction qui appelle les opérateurs de destruction en fonction des poids w1, w2, w3, w4 et w5
function PoidsDestructionVRP(w1, w2, w3, w4, w5, Planning, request_banque, q, pworst, travelTime, TWB, TWE, shiftB, shiftE, meanTaskduration, external_cost, hire_cost, travel_cost, skills, taskByPatient, nurse_workload, P, nurse_by_patient)
    # Normalisation au cas où les poids ne font pas exactement 1
    total = w1 + w2 + w3 + w4 + w5
    w1, w2, w3, w4 = w1 / total, w2 / total, w3 / total, w4/total, w5/total

    r = rand()
    if r < w1
        #println("LE RANDOM A ETE UTILISE")
        #println("RANDOM REMOVE")
        return randomRemove(Planning, request_banque, q, travelTime, TWB, TWE, meanTaskduration, shiftB, shiftE, nurse_workload, nurse_by_patient, taskByPatient)
    elseif r < w1 + w2
        #println("WORST REMOVE")
        return worstRemoval(Planning, request_banque, q, travelTime, TWB, TWE, meanTaskduration, shiftB, shiftE,
                           external_cost, hire_cost, travel_cost, pworst, taskByPatient, nurse_workload, nurse_by_patient)
    elseif r < w1 + w2 + w3
        #println("SHAW REMOVE")
        return relatedRemoval(Planning, request_banque, q, TWB, TWE, travelTime, skills,
        pworst, meanTaskduration, shiftB, shiftE, nurse_workload, nurse_by_patient, taskByPatient)
    elseif r < w1 + w2 + w3 + w4
        #println("WAITING REMOVE")
        return waitingTimeRemoval(Planning, request_banque, travelTime, TWB, TWE, meanTaskduration, shiftB, shiftE, q, pworst, nurse_workload, nurse_by_patient, taskByPatient)
    else
        #println("Concentrated REMOVE")
        return concentratedRemoval(Planning, request_banque, 2, travelTime, TWB, TWE, meanTaskduration, shiftB, shiftE, taskByPatient, nurse_workload, P, nurse_by_patient)
    end
end


# Fonction qui appelle les opérateurs de construction en fonction des poids v1,v2 et v3
function PoidsConstructionVRP(v1, v2, v3, Planning, request_banque, travelTime, TWB, TWE, 
    shiftB, shiftE, meanTaskduration, external_cost, hire_cost, travel_cost, taskByPatient, maxO, skills, nurse_workload, max_weekly_work, rest_hour, nurse_by_patient, maxNumberNurse)

    # Normalisation des poids
    total = v1 + v2 + v3
    v1, v2, v3 = v1 / total, v2 / total, v3 / total

    ran = rand()
    if ran < v1
        return  GreedyInsertion(Planning, request_banque, travelTime, TWB, TWE,
        meanTaskduration, shiftB, shiftE, external_cost, hire_cost, travel_cost,
        taskByPatient, maxO, skills, nurse_workload, max_weekly_work, rest_hour, nurse_by_patient, maxNumberNurse)
    elseif ran < v1 + v2
        return RegretInsertion(Planning, request_banque, travelTime, TWB, TWE,
        meanTaskduration, shiftB, shiftE, external_cost, hire_cost, travel_cost,
        taskByPatient, maxO, skills, nurse_workload, max_weekly_work, 2, rest_hour, nurse_by_patient, maxNumberNurse)
    else
        return RegretInsertion(Planning, request_banque, travelTime, TWB, TWE,
        meanTaskduration, shiftB, shiftE, external_cost, hire_cost, travel_cost,
        taskByPatient, maxO, skills, nurse_workload, max_weekly_work, 3, rest_hour, nurse_by_patient, maxNumberNurse)
    end
end

# Fonction qui applique le critère d'acceptation des solutions dans le LNS
function criteriaVRP(cost_best, cost_new, threshold)
    return cost_new <= threshold*cost_best
end


# Fonction qui retourne le dicotionnaire du nombre d'heure de travail par semaine pour chaque infirmières
function WorkingHourShift(Planning, shiftB, shiftE)
    nurse_hours = Dict{Int, Float64}()

    for day in 1:length(Planning)
        for tournee in Planning[day]
            nurse = tournee[1]
            shift = tournee[2]

            if ismissing(shift)
                duration = 0.0
                if haskey(nurse_hours, nurse)
                    nurse_hours[nurse] += duration
                else
                    nurse_hours[nurse] = duration
                end
                continue  # tournée vide
            end
    
            duration = shiftE[shift] - shiftB[shift]

            if haskey(nurse_hours, nurse)
                nurse_hours[nurse] += duration
            else
                nurse_hours[nurse] = duration
            end
        end
    end

    return nurse_hours
end


# Fontion qui renvoie un dictionnaire avec clé : patient, valeur : Dico{clé : patient, valeur : nombre de fois soigné}
function get_nurses_per_patient(solution, taskbypatient)
    # Dictionnaire : patient_id => Dict(nurse_id => count)
    nurse_counts_by_patient = Dict{Int, Dict{Int, Int}}()

    # Initialiser le dictionnaire pour chaque patient
    for patient_id in keys(taskbypatient)
        nurse_counts_by_patient[patient_id] = Dict{Int, Int}()
    end

    # Parcourir chaque jour de la solution
    for day_routes in solution
        for route_info in day_routes
            if length(route_info) < 3 || isa(route_info[3], Missing)
                continue
            end

            nurse_id = route_info[1]
            path_tasks = route_info[3]
            actual_tasks = [t for t in path_tasks if t != 0]

            for task_id in actual_tasks
                for (patient_id, task_list) in taskbypatient
                    if task_id in task_list
                        if haskey(nurse_counts_by_patient[patient_id], nurse_id)
                            nurse_counts_by_patient[patient_id][nurse_id] += 1
                        else
                            nurse_counts_by_patient[patient_id][nurse_id] = 1
                        end
                        break
                    end
                end
            end
        end
    end

    return nurse_counts_by_patient
end



# Paramètre du LNS :


iteration = 10000 # nombre d'itération du LNS
pworst = 5 # Paramètre rajoutant de l'aléatoire
threshold = 1.05 # critère d'acceptance de la solution
maxO = 0 # maximum d'overtime autorisé


# poids destruction
w1 = 5
w2 = 2
w3 = 1
w4 = 3
w5 = 3

# Poids construction
v1 = 5
v2 = 2
v3 = 1



function LNS_VRP_week(travelTime, TWB, TWE, meanTaskduration, shiftB, 
    shiftE, external_cost, hire_cost, overtime_cost, travel_cost, taskByPatient, maxO, skills, maxNumberNurse, max_weekly_work, rest_hour, 
    iteration, pworst, threshold, N, D, J, lastTasksOnDay, w1, w2, w3, w4, w5, v1, v2, v3)

    # Initialisation de la solution
    Planning, request_banque = InitialSolution(N, D, J, lastTasksOnDay)
    S = Planning

    
    # Dico des heures hebdo des infirmières
    nurse_workload = WorkingHourShift(S, shiftB, shiftE)

    #println(nurse_workload)

    # Dico du nombre d'infirmières par patient
    nurse_by_patient = get_nurses_per_patient(S, taskByPatient)


    # On définit les meilleurs valeur pour S et la request banque
    S_best = S
    best_request_banque = request_banque
    best_nurse_workload = nurse_workload
    best_nurse_by_patient = nurse_by_patient

 
    # Cout de la meilleur solution
    best_cost = compute_cost(S_best, best_request_banque, travelTime, external_cost, hire_cost, overtime_cost, travel_cost, shiftB)


    #println("Solution initiale : ", S_best)
    #println("Request banque initiale : ", best_request_banque)
    #println()

    i = iteration 
    while i > 0
        destruction = rand(10:25)
        q = floor(Int, J * destruction * 0.01)

        #println("valeur de l'itération : ", i)
        S2 = deepcopy(S)
        request_banque2 = deepcopy(request_banque)
        nurse_workload2 = deepcopy(nurse_workload)
        nurse_by_patient2 = deepcopy(nurse_by_patient)

        
        #println("Solution avant destruction : ", S2)
        #println("Request banque de S2 avant destruction : ", request_banque2)
        #println()
        
        
        # On applique un opérateur de destruction 
        PoidsDestructionVRP(w1, w2, w3, w4, w5, S2, request_banque2, q, pworst, travelTime, TWB, TWE, shiftB, shiftE,
        meanTaskduration, external_cost, hire_cost, travel_cost, skills, taskByPatient, nurse_workload2, P, nurse_by_patient2)
        
        #println("Destruction a été faite")

        #println("Solution avant après destruction : ", S2)
        #println("Request banque de S2 après destruction : ", request_banque2)
        #println()

        PoidsConstructionVRP(v1, v2, v3, S2, request_banque2, travelTime, TWB, TWE, shiftB, shiftE, 
        meanTaskduration, external_cost, hire_cost, travel_cost, taskByPatient, maxO, skills, nurse_workload2, 
        max_weekly_work, rest_hour, nurse_by_patient2, maxNumberNurse)

        #println("Reconstruction a été faite")
  
        #println("Solution avant après reconstruction : ", S2)
        #println("Request banque de S2 après reconstruction: ", request_banque2)
        #println()

        current_cost = compute_cost(S2, request_banque2, travelTime, external_cost,
                                   hire_cost, overtime_cost, travel_cost, shiftB)

        if current_cost < best_cost
            S_best = S2
            best_cost = current_cost
            best_request_banque = request_banque2
            best_nurse_workload = nurse_workload2
            best_nurse_by_patient = nurse_by_patient2
        end

        if criteriaVRP(best_cost, current_cost, threshold)
            S = S2
            request_banque = request_banque2
            nurse_workload = nurse_workload2
            nurse_by_patient = nurse_by_patient2
        end

        i = i-1
    end

    return S_best, best_request_banque, best_cost, best_nurse_workload, best_nurse_by_patient
end



#=
using Profile
using FlameGraphs
using ProfileSVG


Planning, request_banque = InitialSolution(N,D,J,lastTasksOnDay)


# Configurer le profilage (plus fréquent)
Profile.init(delay=0.0001)

# Nettoyer les anciens profils
Profile.clear()

# Profiler
@profile LNS_VRP_week(travelTime, TWBeg, TWEnd, meanTaskDuration, shiftBeg, 
shiftEnd, c_e, c_a, c_h, c_o, taskByPatient, maxO, Skills, maxNursesPerPatient, max_weekly_work, rest_hour, iteration, pworst,
threshold, N, D, J, lastTasksOnDay, w1, w2, w3, w4, w5, v1, v2, v3)

Profile.print()





# Générer le flamegraph
fg = flamegraph()

# Sauvegarder si fg est valide
if fg !== nothing
    ProfileSVG.save("flamegraph.svg", fg)
    println("✅ Flamegraph saved to flamegraph.svg")
else
    println("❌ Aucun échantillon collecté, essaye d’augmenter la charge ou ajuster le délai.")
end
=#


#=
@btime begin
    S, B, C, nurse_workload, nurse_by_patient = LNS_VRP_week($travelTime, $TWBeg, $TWEnd, $meanTaskDuration, $shiftBeg, 
    $shiftEnd, $c_e, $c_h, $c_o, $c_t, $tasktopatient, $maxO, $Skills, $maxNursesPerPatient, $max_weekly_work, $rest_hour, $iteration, $pworst,
    $threshold, $N, $D, $J, $lastTasksOnDay, $w1, $w2, $w3, $w4, $w5, $v1, $v2, $v3)
end
=#

tasktopatient = deepcopy(taskByPatient) # j'utilise une copie de la variable taskByPatient parce que dans le LNS je suis amené à modifier cette variable
                                        # Cependant j'ai essayé sans faire de copie et en utilisant taskByPatient et ça marche aussi


S, B, C, nurse_workload, nurse_by_patient = LNS_VRP_week(travelTime, TWBeg, TWEnd, meanTaskDuration, shiftBeg, 
shiftEnd, c_e, c_h, c_o, c_t, tasktopatient, maxO, Skills, maxNursesPerPatient, maxWW*60, R, iteration, pworst,
threshold, N, D, J, lastTasksOnDay, w1, w2, w3, w4, w5, v1, v2, v3)


println("Solution du VRP : ", S)
println()
println("Request_banque du VRP : ", B)
println()

println()
cost, travel, extern, hire, overtime = compute_costbis(S, B, travelTime, c_e, c_h, c_o, c_t, taskByPatient, shiftBeg)
println("Cout total du VRP : ", cost)
println("Cout de travel du VRP : ", travel)
println("Cout d'externalisation du VRP : ", extern)
println("Cout d'hiring du VRP : ", hire)
println("Cout d'overtime du VRP : ", overtime)
println()


println("Heure de travail hebdomadaire des infirmières : ", nurse_workload)
println()
println("Le nombre d'infirmière par patient ", nurse_by_patient)
println()





# Fonction qui retourne un dictionnaire indiquant le nombre d'heure supplémentaire de chaque infirmière dans la solution
function overtime_by_nurse(solution, shiftE)
    overtime = Dict{Int, Float64}()

    for day in solution
        for tour in day
            nurse_id = tour[1]
            shift_id = tour[2]
            G = tour[4]
            if !ismissing(G)
                end_time = G[1][end]
            else
                continue
            end

            # Fin théorique du shift
            if !ismissing(shift_id)
                shift_end = shiftE[shift_id]
            
                # Heures sup
                extra = max(0.0, end_time - shift_end)

                # Ajout au cumul
                if haskey(overtime, nurse_id)
                    overtime[nurse_id] += extra
                else
                    overtime[nurse_id] = extra
                end
            end
        end
    end

    return overtime
end
println("Le nombre d'heures sup pour chaque infirmières dans la semaine : ", overtime_by_nurse(S, shiftEnd))
println()



# Fonction qui renvoie le nombre de tâches effectuées par chaque infirmière dans la solution
function tasks_by_nurse(solution::Vector)
    task_count = Dict{Int, Int}()

    for day in solution
        for tour in day
            nurse_id = tour[1]
            path = tour[3]
            if !ismissing(path)

                # Supposer que le chemin commence et se termine par le dépôt (0), donc on enlève ça
                nb_tasks = count(x -> x != 0, path)

                if haskey(task_count, nurse_id)
                    task_count[nurse_id] += nb_tasks
                else
                    task_count[nurse_id] = nb_tasks
                end
            end
        end
    end

    return task_count
end
println("Le nombre de tache par infirmière par semaine : ", tasks_by_nurse(S))



# Fonction qui affiche la solution de manière plus visible
function affichage(planning, request_banque)
    nb_days = length(planning)
    nb_nurses = 0

    # Dictionnaires pour stocker les shifts et trajets
    nurse_shifts = Dict{Int, Vector{Int}}()
    nurse_routes = Dict{Int, Vector{Vector{Int}}}()

    assigned_tasks_per_day = [Set{Int}() for _ in 1:nb_days]

    for (day_idx, day) in enumerate(planning)
        for entry in day
            if entry === nothing
                continue
            end

            nurse_id, shift, route = entry

            route = route === nothing ? Int[] : route
            nb_nurses = max(nb_nurses, nurse_id)

            if !haskey(nurse_shifts, nurse_id)
                nurse_shifts[nurse_id] = fill(0, nb_days)
                nurse_routes[nurse_id] = [Int[] for _ in 1:nb_days]
            end

            nurse_shifts[nurse_id][day_idx] = shift === missing ? 0 : Int(shift)


            if route === nothing || route === missing
                route = Int[]
            end
            
            route_cleaned = length(route) ≥ 2 ? route[2:end-1] : Int[]
            
            nurse_routes[nurse_id][day_idx] = route_cleaned
            union!(assigned_tasks_per_day[day_idx], route_cleaned)
        end
    end

    

    # Extraction de toutes les tâches demandées, jour par jour
    all_tasks_per_day = [Set{Int}() for _ in 1:nb_days]
    for (day_idx, day_requests) in enumerate(request_banque)
        if day_requests !== nothing && day_requests isa AbstractVector
            union!(all_tasks_per_day[day_idx], day_requests)
        end
    end

    #println(all_tasks_per_day)


    # Calcul des patients externalisés
    outsourced_patients_per_day = Vector{Vector{Int}}()
    for day_idx in 1:nb_days
        externalized = sort(collect(setdiff(all_tasks_per_day[day_idx], assigned_tasks_per_day[day_idx])))
        push!(outsourced_patients_per_day, externalized)
    end


    # Affichage
    print("Outsourced tasks : ")
    for day_ext in outsourced_patients_per_day
        print("| ", join(day_ext, " "), " ")
    end
    println("\n")

    for nurse_id in 1:nb_nurses
        shifts = get(nurse_shifts, nurse_id, fill(0, nb_days))
        routes = get(nurse_routes, nurse_id, [Int[] for _ in 1:nb_days])
        print("Nurse $nurse_id : ", shifts, " | ")
        for day_route in routes
            print(day_route, ", ")
        end
        println()
    end
end

println()
println("AFFICHAGE DE LA SOLUTION")


println()
println(affichage(S, B))
println()




# Checker de solution. Il check notamment si : 
    # Chaque infirmière possède la compétence pour effectuer la tâche
    # Il n'y pas de doublon dans la solution
    # Les fenêtres horaires sont bien respectées
    # Les chemins partent et arrivent bien au dépôt
    # Les infirmières ne dépassent pas le temps de travail hebdomadaire
    # Les infirmières ont un temps de repos suppérieur au temps de repos minimum 
    # Toutes les tâches sont présentes dans la solution
    # La continuité des soins est bien respecetée

function CheckSolution(Planning, request_banque, travelTime, TWB, TWE, shiftB, shiftE, meanTaskduration, Skills, TaskbyPatient, maxWW, rest_hour, maxO, maxNursesPP, N, J)
    workload = zeros(N)
    shift_info = zeros(N, length(Planning))
    Task = deepcopy(request_banque)

    # --- Dico Suivi : patient → (nurse → nb de soins) ---
    Suivi = Dict{Int, Dict{Int, Int}}()

    println()
    jour = 1
    for day in Planning
        for route in day
            nurse = route[1]
            shift = route[2]
            path = route[3]
            
            if ismissing(shift)
                continue
            end

            shift_info[nurse, jour] = shift
            nurse_skill = Skills[nurse]

            for k in 2:length(path)-1
                current_task = path[k]

                # --- Check des compétences ---
                if !(current_task in nurse_skill)
                    println("La contrainte des skills n'est pas respectée pour la tâche ", current_task, " et l'infirmière ", nurse)
                    return false
                end

                # --- Check des doublons dans Task ---
                if !(current_task in Task)
                    push!(Task, current_task)
                else
                    println("La tâche ", current_task, " apparait deux fois dans la solution : solution erronée")
                    return false
                end

                # --- Remplissage du dico Suivi ---
                # Trouver le patient correspondant à cette tâche
                for (patient, t_list) in TaskbyPatient
                    if current_task in t_list
                        if !haskey(Suivi, patient)
                            Suivi[patient] = Dict{Int, Int}()
                        end
                        Suivi[patient][nurse] = get(Suivi[patient], nurse, 0) + 1
                        break
                    end
                end
            end

            # --- Gestion du temps de travail ---
            duration = shiftE[shift] - shiftB[shift]
            workload[nurse] += duration

            # --- Check des fenêtres horaires ---
            t = shiftB[shift]
            for k in 2:length(path)-1
                current_task = path[k]
                t = max(TWB[current_task], t + travelTime[path[k-1]+1, current_task+1])
                if t > TWE[current_task]
                    println("Fenêtre horaire non respectée pour la tâche ", current_task,
                            " arrivée à ", t, " alors que la fenêtre finit à ", TWE[current_task])
                    return false
                end
                t += meanTaskduration[current_task]
            end

            # --- Vérification retour dépôt ---
            t += travelTime[path[end-1]+1, path[end]+1]
            if t > shiftE[shift] + maxO
                println("Le shift est dépassé pour la nurse ", nurse, " retour à ", t, " au lieu de ", shiftE[shift])
                return false
            end
        end
        jour += 1
    end

    println("Les contraintes de compétences sont bien respectées")
    println()
    println("Aucune tâche n'apparait deux fois dans la solution")
    println()
    println("Les contraintes de fenêtres horaires sont bien respectées\n")

    # --- Check temps de travail hebdo ---
    for k in 1:length(workload)
        if workload[k] > maxWW
            println("La contrainte du temps de travail hebdo n'est pas respectée par la nurse ", k, " qui travaille ", workload[k])
            return false
        end
    end
    println("Les contraintes des temps de travail sont respectées : workload = ", workload, "\n")

    # --- Check temps de repos ---
    for i in 1:size(shift_info, 1)
        for j in 1:size(shift_info, 2)-1
            if shift_info[i,j] > 0 && shift_info[i,j+1] > 0
                rest_time = (24*60 - shiftE[Int(shift_info[i,j])]) + shiftB[Int(shift_info[i,j+1])]
                if rest_time < rest_hour*60
                    println("La contrainte des temps de repos n'est pas respectée pour la nurse ", i, " (repos = ", rest_time, " min)")
                    return false
                end
            end
        end
    end
    println("Les contraintes de temps de repos sont respectées\n")

    # --- Check si toutes les tâches sont présentes ---
    total = 0
    for el in Task
        if el isa AbstractVector
            total += length(el)
        else
            total += 1
        end
    end
    
    if total !== J
        println("Le nombre de tâches présentes est ", total, " au lieu de ", J)
        return false
    end
    println("Toutes les tâches sont présentes : solution valide\n")

    # --- Check maxNursesPP ---
    for (patient, nurse_counts) in Suivi
        if length(keys(nurse_counts)) > maxNursesPP[patient]
            println("Le patient ", patient, " est suivi par ", length(keys(nurse_counts)),
                    " infirmières différentes (max autorisé = ", maxNursesPP, ")")
            return false
        end
    end
    println("Contrainte maxNursesPP respectée\n")

    return true
end


println("CHECK DE LA SOLUTION")
check = CheckSolution(S, B, travelTime, TWBeg, TWEnd, shiftBeg, shiftEnd, meanTaskDuration, Skills, taskByPatient, maxWW*60, R, maxO, maxNursesPerPatient, N, J)
println()
println("la valeur du check : ", check)






#=
function municipalitiesMatrix(travelTime, taskMunicipality)
    nbMun = maximum(taskMunicipality) + 1
    munMatrix = fill(-1.0, nbMun, nbMun)

    for i in 1:length(taskMunicipality)
        for j in i+1:length(taskMunicipality)
            m1 = taskMunicipality[i] + 1  # décalage car les municipalités vont de 0 à N-1
            m2 = taskMunicipality[j] + 1

            # Déterminer les indices dans travelTime avec le décalage dépôt
            idx1 = i + 1
            idx2 = j + 1

            if munMatrix[m1, m2] == -1
                munMatrix[m1, m2] = travelTime[idx1, idx2]
                munMatrix[m2, m1] = travelTime[idx1, idx2]
            end
        end
    end

    for i in 1:nbMun
        munMatrix[i, i] = 3.0
    end

    return munMatrix
end

println()
munTravelTime = municipalitiesMatrix(travelTime, taskMunicipality)
println("La matrice de distance entre municipalitées : ", munTravelTime)




function distance_matrix_to_coordinates(dist_matrix::Matrix{Float64})
    size(dist_matrix, 1) == size(dist_matrix, 2) || error("La matrice doit être carrée")
    
    n = size(dist_matrix, 1)
    H = I - fill(1/n, n, n)
    D_squared = dist_matrix.^2
    B = -0.5 * H * D_squared * H
    
    eigenvalues, eigenvectors = eigen(Symmetric(B))
    idx = sortperm(eigenvalues, rev=true)[1:2]
    
    coords = eigenvectors[:, idx] * Diagonal(sqrt.(eigenvalues[idx]))
    
    return coords
end




coordinates = distance_matrix_to_coordinates(munTravelTime)

println()
println("Coordonnées des points :")
for (i, coord) in enumerate(eachrow(coordinates))
    println("Point $i : (x=$(round(coord[1], digits=2)), y=$(round(coord[2], digits=2)))")
end

println()

# Calcul du centroïde
depot_x = mean(coordinates[:, 1])
depot_y = mean(coordinates[:, 2])
println("Coordonnées estimées du dépôt : (x=$(round(depot_x, digits=2)), y=$(round(depot_y, digits=2)))")





function affichage_graphique(path, coords, depot_x, depot_y, taskMunicipalities)

    # Créer le graphique
    scatter(title = "Chemin parcouru", xlabel = "x", ylabel = "y", legend = false)

    # Tracer les municipalités (cercles colorés par municipalité)
    nbTasks = size(coords, 1)
    for i in 1:nbTasks
        mun_id = taskMunicipalities[i] + 1  # +1 si municipalités sont indexées de 0
        scatter!([coords[i, 1]], [coords[i, 2]], label = false, markersize = 8, color = palette(:tab10)[mun_id])
        annotate!(coords[i, 1], coords[i, 2], text(string(i), :left, 8))
    end

    # Tracer le dépôt
    scatter!([depot_x], [depot_y], markersize = 10, color = :black, marker=:star5, label = "Dépôt")
    annotate!(depot_x, depot_y, text("Dépôt", :bottom, 10))

    # Construire le chemin en coordonnées
    x_coords = [depot_x]  # départ depuis le dépôt
    y_coords = [depot_y]
    for i in path
        mun = taskMunicipalities[i] + 1
        #println(mun)
        push!(x_coords, coords[mun, 1])
        push!(y_coords, coords[mun, 2])
    end
    push!(x_coords, depot_x)  # retour au dépôt
    push!(y_coords, depot_y)

    # Tracer le chemin
    plot!(x_coords, y_coords, lw=2, color=:red, label = false, arrow=true)

    display(current())

    gui()               # force l'ouverture de la fenêtre de graphique

    return
end


function affichage_nurse_path(Planning, nurse, coords, depot_x, depot_y, taskMunicipalities)
    for day in Planning
        path = day[nurse][3]
        reel_path = path[2:end-1]
        #println(reel_path)
        affichage_graphique(reel_path, coords, depot_x, depot_y, taskMunicipalities)
    end
end
        

#path = [78, 90, 94, 93, 91, 81]
#affichage_graphique(path, coordinates, depot_x, depot_y, taskMunicipality)
#affichage_nurse_path(S, 2, coordinates, depot_x, depot_y, taskMunicipality)

=#

