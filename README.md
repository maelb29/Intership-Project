# 🏥 Intership-Project  
**LNS metaheuristic for VRPTW (Vehicle Routing Problem with Time Windows)**

Ce projet a été réalisé dans le cadre d’un **stage de 2e année**.  
Il a pour objectif d’**optimiser des tournées de soins à domicile** en utilisant une métaheuristique de type **Large Neighborhood Search (LNS)**, testée sur les instances proposées par **Paola Cappanera**.

---

## 📂 Organisation des fichiers

- **Implementation VRP best.jl** : Fichier **principal** qui exécute le LNS et retourne la solution.  
- **Operateurs VRP best.jl** : Contient les opérateurs de destruction et de réparation.  
- **Variables Globales VRP best.jl** : Définit les variables et paramètres globaux utilisés par le LNS.  
- **instances/** : Répertoire contenant l’ensemble des **instances de Paola Cappanera**.

---

## 🚀 Comment exécuter le projet ?

### 1. Inclure les fichiers nécessaires

Dans `Implementation VRP best.jl`, ajoutez les lignes suivantes (en adaptant les chemins si nécessaire) :  

```julia
include("Operateurs VRP best.jl")
include("Variables Globales VRP best.jl")
include("instances/nom_instance.jl")



# --------------------------
# Paramètres principaux LNS
# --------------------------

Une fois les chemins des fichiers renseigner on peut changer si on le souhaite les paramètres du LNS qui sont déjà renseignés par défaut. On peut notamment jouer sur : le nombre d'iteration, le paramètre pworst qui applique un biais aléatoire utile pour les opérateurs de destruction, le paramètre threshold qui intervient dans le critère d'acceptation d'une nouvelle solution, le paramètre maxO qui indique le nombre d'heures supplémentaires autorisées, et enfin les paramètres concernant les poids attribués à chaque opérateurs. 


iteration = 10000        # Nombre d'itérations du LNS
pworst = 5               # Biais aléatoire (utilisé par les opérateurs de destruction)
threshold = 1.05         # Critère d'acceptation d'une nouvelle solution
maxO = 0                 # Maximum d'heures supplémentaires autorisées

# --------------------------
# Poids des opérateurs
# --------------------------

# Poids des opérateurs de destruction
w1 = 5
w2 = 2
w3 = 1
w4 = 3
w5 = 3

# Poids des opérateurs de construction
v1 = 5
v2 = 2
v3 = 1

Une fois tout cela effectué on peut lancer le fichier Implementation VRP best.jl et celui retournera la solution du problème.