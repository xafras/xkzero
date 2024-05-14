# -*- coding: utf-8 -*-
#   Xafras xyZw Krazos
#   Version 1.0
#   Ouvert le : Sat Jan 14 02:20:38 2023
#   Fermé le : Sat Jan 14 02:20:38 2023

#   Initialisation

""" Ce package implemente la methode de Newton du premier ordre et la recherche
dichotomique de zero pour une fonction reelle dont l'espace d'arrive est muni
d'une mesure quelconque. Cette methode ne prend pas en charge les fonctions a
variables non reelles. """

#   Fonction

def Newton(x0, get_f, get_fp, get_foverfp, get_mu=abs, eps=10**(-3), N=1000, renvoyer_liste=False, debug=False):
    """ Methode de Newton pour approcher un zero proche de x0 de la fonction f
    de derivee fp et telle que foverfp = f/f'. Le zero est approche a eps pres
    selon la mesure mu.
    
    Parametres
    ----------
    x0 : float
        Valeur initiale, une valeur proche du resultat attendu permet un
        meilleur resulat.
    get_f : 'a -> 'b 
        Fonction f a annuler.
    get_fp : 'a -> 'b
        Derivee de la fonction f a annuler.
    get_foverfp : 'a -> 'b
        Expression de f/f' simplifiee.
    get_mu : 'a -> float
        Mesure sur l'image de f. Valeur absolue par defaut.
    eps : float, optionnel
        Tolerance du defaut de nullite de f(x). 10**(-3) par defaut.
    N : int, optionnel
        Limite maximale d'operation avant abandon. 1000 par defaut.
    renvoyer_liste : bool, optionnel
        Quand definit a True, la methode renvoie la liste de toute les valeurs
        prise par la suite. False par defaut.
    debug : bool, optionnel
        Quand definit a True, les informations de chaque etapes de calcul sont
        affichees en direct dans la console. False par defaut.

    Resultats
    -------
    x : float
        Zero de f approche a eps pres.
    mu(f(x)) : float
        Defaut de nullite de f(x).
    n!=N : bool
        Booleen valant True si le calcul s'est acheve correctement et False
        s'il a ete abrege.
    n : int
        Nombre d'etapes de calcul realisees.
    L : float list
        Zeros de f a chaque etape de l'agorithme. """
    n = 0
    x = x0
    L = [x]
    while n < N and get_mu(get_f(x)) > eps and get_fp(L[-1]) != 0:
        n += 1
        x = L[-1]
        L.append(x - get_foverfp(x))
        if debug:
            print(f'Newton: n={n}, x={x}')
    if renvoyer_liste:
        return L,get_mu(get_f(x)),n!=N,n
    return x,get_mu(get_f(x)),n!=N,n

def dichotomie(xmin, xmax, get_f, get_mu=abs, eps=10**(-3), N=100, renvoyer_liste=False, debug=False):
    """ Methode dichotomique de recherche de zero de la fonction f compris
    entre xmin et xmax. Le zero est approche a eps pres selon la mesure mu.

    Parameters
    ----------
    xmin : float
        DESCRIPTION.
    xmax : float
        DESCRIPTION.
    get_f : float -> 'a
        DESCRIPTION.
    get_mu : 'a -> float
        Mesure sur l'image de f. Valeur absolue par defaut.
    eps : float, optionnel
        Tolerance du defaut de nullite de f(x). 10**(-3) par defaut.
    N : int, optionnel
        Limite maximale d'operation avant abandon. 1000 par defaut.
    renvoyer_liste : bool, optionnel
        Quand definit a True, la methode renvoie la liste de toute les valeurs
        prise par la suite. False par defaut.
    debug : bool, optionnel
        Quand definit a True, les informations de chaque etapes de calcul sont
        affichees en direct dans la console. False par defaut.

    Resultats
    -------
    x : float
        Zero de f approche a eps pres.
    mu(f(x)) : float
        Defaut de nullite de f(x).
    n!=N : bool
        Booleen valant True si le calcul s'est acheve correctement et False
        s'il a ete abrege.
    n : int
        Nombre d'etapes de calcul realisees.
    L : float list
        Zeros de f a chaque etape de l'agorithme. """
    n = 0
    m = (xmax+xmin)/2
    fm = get_f(m)
    L = [m]
    while n < N and get_mu(fm) > eps :
        n += 1
        if fm == 0 :
            break
        elif fm < 0 :
            xmin = m
        else :
            xmax = m
        m = (xmax+xmin)/2
        fm = get_f(m)
        L.append(m)
        if debug:
            print(f'Dichotomie: n={n}, x={m}')
    if renvoyer_liste:
        return L,get_mu(get_f(L[-1])),n!=N,n
    return L[-1],get_mu(get_f(L[-1])),n!=N,n

#   Exemple

if __name__=="__main__":    #   exemple avec x -> x**2 - 1 (suite de Heron)
    
    #   Newton
    def get_f(x):
        return x**2 - 2
    
    def get_fp(x):
        return 2*x
    
    def get_foverfp(x):
        if x==0:
            raise "Erreur : Division par zero :c"
        return x/2 - 1/x
    
    def get_mu(x):
        return abs(x)
    
    print(Newton(1, get_f, get_fp, get_foverfp, get_mu,debug=(True)))
    
    #   Dichotomie
    print(dichotomie(1, 4, get_f, debug=(True),renvoyer_liste=(True)))

#   Fin