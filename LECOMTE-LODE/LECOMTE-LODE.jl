# Soline LECOMTE - Gaël LODE

using DelimitedFiles, Statistics# On va utiliser JuMP et GLPK


# Structures contenant les données d'un problème de bin-packing mono-dimensionnel
struct typeData
	m::Int64 # Nombre de machine
	n::Int64 # Nombre de tache
	d::Vector{Int64} #Durée de chacune des taches
end

struct typeLSA
	maxLSA::Int64
	borneMax::Float64
	borneMoy::Float64
	ratio::Float64
end
struct typeLPT
	maxLPT::Int64
	ratio::Float64
end
struct typeMyAlgo
	maxMA::Int64
	ratio::Float64
end
# Fonction prenant en entrant un nom de fichier et retournant les données du problème
function parserFichier(nomFichier::String)
    # Ouverture d'un fichier en lecture
    f = open(nomFichier,"r")    
    s= readdlm(f,':',Int64,'\n')
    inm = deepcopy(s[1])
    inn = deepcopy(s[2])
    tab::Vector{Int64} = Vector{Int64}(undef,inn)
    for i in 1:inn
    	tab[i] = deepcopy(s[i+2])
    end
    return typeData(inm, inn, tab)
end


function parserClavier(clavier::String)
	i::Int64 = 1
	while i <= length(clavier) && clavier[i] != ':'
		i = i+1
	end
    inm = parse(Int64,deepcopy(clavier[1:i-1]))
   
    i=i+1
    tmp1 = i
    while i <= length(clavier) && clavier[i] != ':'
    	i=i+1
   	end
   	inn = parse(Int64,deepcopy(clavier[tmp1:i-1]))
   	i=i+1
   	tmp2 = i
   	tab::Vector{Int64} = Vector{Int64}(undef,inn)
   	j= 1
   	while i <= length(clavier) && clavier[i] != '\n'&& j <= inn
   			if(clavier[i] == ':')
   				tab[j] = parse(Int64,deepcopy(clavier[tmp2:i-1]))
   				tmp2 = i+1
   				j=j+1
   			end
   			i = i+1
   	end
   	tab[j] = parse(Int64,deepcopy(clavier[tmp2:i-1]))
    return typeData(inm,inn,tab)
end

function parser1Instance(m::String)
	inm = parse(Int64,m)
	inn = inm*2 +1 
	i=1
	tab::Vector{Int64} = Vector{Int64}(undef,inn)
	for i in 1:3
		tab[i] = inm
	end
	i=i+1
	j = 1
	while i <= inn
		tab[i] = inm +j
		tab[i+1] = inm+j
		i=i+2
		j=j+1
	end
	return typeData(inm,inn,tab)
end

function parserkInsances(mm::Int64, nn::Int64, k::Int64, mini::Int64, maxi::Int64)
	sortie::Vector{typeData} = Vector{typeData}(undef,k)
	for i in 1:k
		inm = mm
		inn = nn
		tab::Vector{Int64} = Vector{Int64}(undef,inn)
		for j in 1:nn
			tab[j] = round(rand(mini:maxi))
		end
		sortie[i] = typeData(inm,inn,tab)
	end
	return sortie
end


function passageAlgoSimple(donne::typeData)
#LSA
MLSA = fill(0,donne.m)

for i in 1:donne.n
	mini = 1
	for j in 1:donne.m
		if  MLSA[j] < MLSA[mini]
			mini = j
		end
	end
	MLSA[mini] += donne.d[i]
end

resLSA = maximum(MLSA)
borneInfMax = maximum(donne.d)
borneInfMoy = sum(donne.d)/donne.m
maxb = max(borneInfMax, borneInfMoy)
res1=typeLSA(resLSA, borneInfMax, borneInfMoy, resLSA/maxb)


#LPT
MLPT= fill(0,donne.m)
tab1 = deepcopy(donne.d)
sort!(tab1,rev=true)
for i in 1:donne.n
	mini = 1
	for j in 1:donne.m
		if MLPT[j] < MLPT[mini]
			mini = j
		end
	end
	MLPT[mini] += tab1[i]
end
resLPT = maximum(MLPT)
res2=typeLPT(resLPT, resLPT/maxb)

#My Algo
#On trie les durée dans l'ordre croissant puit on fait la même chose qu'avec le LPT
MMA= fill(0,donne.m)
tab2 = deepcopy(donne.d)
sort!(tab2)
for i in 1:donne.n
	mini = 1
	for j in 1:donne.m
		if MMA[j] < MMA[mini]
			mini = j
		end
	end
	MMA[mini] += tab2[i]
end
resMA = maximum(MMA)
res3=typeMyAlgo(resMA, resMA/maxb)

println("résultat obtenu\n
	Borne inférieur 'maximum' = $(res1.borneMax) \n
	Borne inférieur 'moyenne' = $(res1.borneMoy) \n
	Résultat LSA = $(res1.maxLSA) \n
	ratio LSA = $(res1.ratio)\n
	Résultat LPT = $(res2.maxLPT) \n
	ratio LPT = $(res2.ratio)\n
	Résultat My Algo = $(res3.maxMA) \n
	ratio My Algo = $(res3.ratio)")
end

function passageAlgoK(donne::Vector{typeData}, k::Int64)
#LSA
ratioMLSA=fill(0.0,k)
maxbLSA = fill(0.0,k)
for l in 1:k
	MLSA = fill(0,donne[l].m)
	for i in 1:donne[l].n
		mini = 1
		for j in 1:donne[l].m
			if  MLSA[j] < MLSA[mini]
				mini = j
			end
		end
		MLSA[mini] += donne[l].d[i]
	end

	resLSA = maximum(MLSA)
	borneInfMax = maximum(donne[l].d)
	borneInfMoy = sum(donne[l].d)/donne[l].m
	maxb = max(borneInfMax, borneInfMoy)
	maxbLSA[l] = maxb
	ratioMLSA[l]= resLSA/maxb
end

#LPT
ratioMLPT=fill(0.0,k)
for l in 1:k
	MLPT= fill(0,donne[l].m)
	tab1 = deepcopy(donne[l].d)
	sort!(tab1,rev=true)
	for i in 1:donne[l].n
		mini = 1
		for j in 1:donne[l].m
			if MLPT[j] < MLPT[mini]
				mini = j
			end
		end
		MLPT[mini] += tab1[i]

	end
	resLPT = maximum(MLPT)
	ratioMLPT[l] = resLPT/maxbLSA[l]
end
#My Algo
#On trie les durée dans l'ordre croissant puit on fait la même chose qu'avec le LPT
ratioMMA=fill(0.0,k)
for l in 1:k
	MMA= fill(0,donne[l].m)
	tab2 = deepcopy(donne[l].d)
	sort!(tab2)
	for i in 1:donne[l].n
		mini = 1
		for j in 1:donne[l].m
			if MMA[j] < MMA[mini]
				mini = j
			end
		end
		MMA[mini] += tab2[i]

	end
	resMMA = maximum(MMA)
	ratioMMA[l] = resMMA/maxbLSA[l]
end

println("résultat obtenu\n
	ratio LSA = $(mean(ratioMLSA))\n
	ratio LPT = $(mean(ratioMLPT))\n
	ratio My Algo = $(mean(ratioMMA))")
end

function fichier()
	println("Mode : Depuis un fichier \n
			Le fichier doit etre de la forme 'm:n:d1:...dn'")
	userinp = readline(stdin)
	donne= parserFichier(userinp)
	return passageAlgoSimple(donne)
end

function clavier()
	println("Mode : au clavier \n
			l'entrée doit etre de la forme 'm:n:d1:...dn'")
	userinp = readline(stdin)
	donne::typeData = parserClavier(userinp)
	return passageAlgoSimple(donne)
end

function uneInstance()
	println("Mode : Génération d'une instance \n
		indiquer le nombre de machine m")
	userinp = readline(stdin)
	donne::typeData = parser1Instance(userinp)
	return passageAlgoSimple(donne)
end

function kInstances()
	println("Mode : Génération k instances\n
indiquer le nombre de machine m")
	userinp = readline(stdin)
	m = parse(Int64, userinp)
	println("indiquer le nombre de tache n")
	userinp = readline(stdin)
	n = parse(Int64, userinp)
	println("indiquer le nombre d'instance k")
	userinp = readline(stdin)
	k = parse(Int64, userinp)
	println("indiquer la durée min")
	userinp = readline(stdin)
	mini = parse(Int64, userinp)
	println("indiquer la durée max")
	userinp = readline(stdin)
	maxi = parse(Int64, userinp)
	donne::Vector{typeData} = parserkInsances(m,n,k,mini,maxi)
	return passageAlgoK(donne,k)

end


#fonction principal pour les objet 1D
function main()

	println("Selectionner un mode d'entrée\n
			1 - Depuis un fichier\n
			2 - Au clavier\n
			3 - Génération d'une instance\n
			4 - Génération k instances")
	userinp::String = readline(stdin)
	if userinp == "1" # il n'y a pas de switch avec julia :)
		fichier()
	elseif userinp == "2"
		clavier()
	elseif userinp == "3"
		uneInstance()
	elseif userinp == "4"
		kInstances()
	else 
		println("Entrée incorecte")
		return 1
	end
end