########################################### EXERCICE 1 ###########################################

# QUESTION 1 ----

# PARAMETRE

n <- 100
Sigma <- 0.2


    # a/ SIMULATION DES X_i QUI SUIVENT UNE LOI UNIFORME SUR [0,1] ----

    X <-  runif(n, min = 0, max = 1)
    X

    # b/ SIMULATION DES Y_i QUI SUIVENT UNE LOI UNIFORME SUR [0,1] ----

    # ETAPE 1: ON DEFINIT D'ABORD LA FONCTION f
    
    f <- function(x){
      return((x^2*2^(x-1)-(x-0.5)^3)*sin(10*x))
    }
    
    # ETAPE 2: ON SIMULE ENSUITE LE BRUIT Epsilon
    
    Epsilon = rnorm(n, mean = 0, sd = 1)
    Epsilon
    
    # ETAPE 3: SIMULATION DES Y_i
    
    Y = f(X) + Sigma * Epsilon 
    Y

# QUESTION 2 ----

    # install.packages("ggplot2")
    library(ggplot2)
    
    # CREATION D'UNE DATAFRAME QUI CONTIENT LES DONNEES SIMULEES
    
    Data_X_Y <- data.frame(X = X, Y = Y)
    Data_X_Y
    
    # REPRESENTATION DU NUAGE DE POINT GRACE A LA FONCTION ggplot
    
    ggplot(Data_X_Y, aes(x = X, y = Y)) + 
      
      geom_point(aes(color = "Simulation")) + 
      
      stat_function(fun = f, aes(color = "Fonction f"), size = 1.5) +
      
      labs(
        
        title = "Simulation des données et fonction f",
        
        x = "X",
        
        y = "Y",
        
        color = "Légende"
        
      ) +
      
      theme_minimal() + 
      
      theme(
        
        plot.title = element_text(hjust = 0.5),
        
        legend.position = "top"
        
      )

# QUESTION 3 ----

    # a/ LES 5 PREMIERES FONCTIONS DE LA BASE DE FOURIER ----

    # ON POSE k = 2j
    
    base_trigonometrique <- function(x, k) {
      
      if (k == 1) {
        
        return(1)
        
      }
      
      if (k %% 2 == 0) {
        
        return(sqrt(2) * cos(pi * k * x))
        
      }
      
      return(sqrt(2) * sin(pi * (k - 1) * x))
      
    }
    
    print(base_trigonometrique(5, 1)) # Test
    
    # CREATION D'UNE SEQUENCE DE VALEURS t
    
    t <- seq(from = 0, to = 1, length = 1000)
    
    # CREATIOON D'UNE DATAFRAME QUI CONTIENT LES VALEURS DE t ET LES FONCTION CORRESPONDANTES
    
    data <- data.frame(
      
      t = t,
      
      base_trigonometrique_1 = base_trigonometrique(t, 1),
      
      base_trigonometrique_2 = base_trigonometrique(t, 2),
      
      base_trigonometrique_3 = base_trigonometrique(t, 3),
      
      base_trigonometrique_4 = base_trigonometrique(t, 4),
      
      base_trigonometrique_5 = base_trigonometrique(t, 5)
      
    )
    
    # REPRESENTATION GRAPHIQUE
    
    ggplot(data, aes(x = t)) +
      
      geom_line(aes(y = base_trigonometrique_1), color = "red") +
      
      geom_line(aes(y = base_trigonometrique_2), color = "brown") +
      
      geom_line(aes(y = base_trigonometrique_3), color = "purple") +
      
      geom_line(aes(y = base_trigonometrique_4), color = "green") +
      
      geom_line(aes(y = base_trigonometrique_5), color = "pink") +
      
      ylim(-2, 2) +
      
      labs(
        
        y = "base_trigonometrique",
        title = "Les 5 premières fonctions de la base trigonométrique",
        
        main = "Les 5 premières fonctions de la base trigonométrique"
        
      )

    # b/ CONSTRUCTION DE L'ESTIMATEUR PAR PROJECTION DE f ----

    # PARAMETRE
    
    N = c(5, 10, 15, 20, 30, 50, 70)
    
    # ON DEFINIT D'ABORD  LA FONCTION Y(BASE TRIGONOMETRIQUE)
    
    Y_base_trigonometrique = function(j){
      somme = 0
      for (i in 1:n){
        somme = somme + Y[i]*base_trigonometrique(X[i], j)
      }
      return(somme/n)
    }
    
    # MATRICE DE PROJECTION G 
    
    G = function(j, k){
      somme = 0
      for(i in 1:n){
        somme = somme + base_trigonometrique(X[i], j)*base_trigonometrique(X[i], k)
      }
      return(somme/n)
    }
    
    # L'ESTIMATEUR DE TETA
    
    Teta_chapeau = function(j, N){
      TETA = solve(G(j, j), Y_base_trigonometrique(j))
      return(TETA)
    }
    
    # L'ESTIMATEUR DE f
    
    f_chapeau = function(X, N){
      somme = 0
      for (i in 1:N){
        somme = somme + Teta_chapeau(i, N)*base_trigonometrique(X, i)
      }
      return(somme)
    }
    
    # REPRESENTATION GRAPHIQUE DE L'ESTIMATEUR
    
    t = seq(from = 0, to = 1, length = 1000)
    
    data = data.frame(t = t, f = f(t), f_chapeau = sapply(t, f_chapeau, N = N)) 
    
    ggplot(data, aes(x = t)) +
      
      geom_line(aes(y = f), color = "blue", linetype = "solid") +
      
      geom_line(aes(y = f_chapeau), color = "red", linetype = "dashed") +
      
      labs(title = "Fonction f et f_chapeau", x = "t", y = "Valeur") +
      
      theme_minimal()

    # c/ LA VALEUR DE N OPTIMALE PAR VALIDATION CROISEE ----

    # L'ESTIMATEUR DE f sans X[i]
    
    f_chapeau_moins_i = function(X, N, i){
      somme = 0
      for(j in 1:N){
        somme = somme + Teta_chapeau(j, N)* base_trigonometrique(X, j)*(j != i)
      }
      return(somme)
    }
    
    # LE CRITERE A MINIMISER
    
    critere = function(N){
      somme = 0
      for(i in 1:n){
        somme = somme + (Y[i] -  f_chapeau_moins_i(X[i], N, i))^2
      }
      return(somme)
    }
    
    # CALCUL DE N OPTIMAL ET REPRESENTATION GRAPHIQUE
    
    N_optimal = optimize(critere, 0:50)$minimum
    N_optimal 
    
    t = seq(from = 0, to = 1, length = 1000)
    
    data = data.frame(t = t, f = f(t), f_chapeau = sapply(t, f_chapeau, N = N_optimal)) 
    
    ggplot(data, aes(x = t)) +
      
      geom_line(aes(y = f), color = "blue", linetype = "solid") +
      
      geom_line(aes(y = f_chapeau), color = "red", linetype = "dashed") +
      
      labs(title = "Fonction f et f_chapeau", x = "t", y = "Valeur") +
      
      theme_minimal()

# QUESTION 4 ----
    
    # PARAMETRE 
    
    n = 100
    N_simulation = 200
    sigma = 0.2
    N_result = rep(0, N_simulation)
  
    # ON SIMULE N FOIS
    
    for(i in 1:N_simulation){
      X = runif(n, 0, 1)
      eps = rnorm(n, 0, 1)
      Y = f(X) + sigma*eps
      N_result[i] = optimize(critere, 0:50)$minimum
    }
    hist(N_result, main = 'Histogramme des 200 D optimaux', prob = TRUE)

########################################### EXERCICE 2 ###########################################
    
    # IMPORTATION DE LA BASE DE DONNEES ----
    ozone <- read.table("https://r-stat-sc-donnees.github.io/ozone.txt",header=T)
 
    # VARIABLES CHOISIES ----
    
    X = ozone$T12
    Y = ozone$maxO3

    # ON NORMALISE X
    X = (X - min(X))/(max(X) - min(X))
    
    # REPRESENTATION GRAPHIQUE ----
    # CREATION D'UNE DATAFRAME QUI CONTIENT LES DONNEES SIMULEES
    
    Data_X_Y <- data.frame(X = X, Y = Y)
    Data_X_Y
    

    # REPRESENTATION DU NUAGE DE POINT GRACE A LA FONCTION ggplot

    ggplot(ozone, aes(x = X, y = Y)) + 
      
      geom_point(aes(color = "Donnees d'ozone"))  +
      
      labs(
        
        title = "Nuage de points de Y en fonction de X",
        
        x = "T12",
        
        y = "max03",
        
        color = "Légende"
        
      ) +
      
      theme_minimal() + 
      
      theme(
        
        plot.title = element_text(hjust = 0.5),
        
        legend.position = "top"
        
      )
    
    # HISTOGRAMME DE X
    
    hist(X, main = 'Histogramme de la température', xlab = 'Température')
    

    # PARAMETRE
    N = c(5, 10, 15, 20, 30, 50, 70)
    n = as.numeric(length(X))
    

    # ESTIMATION ----
    
    # ON DEFINIT D'ABORD  LA FONCTION Y(BASE TRIGONOMETRIQUE)
    
    Y_base_trigonometrique = function(j){
      somme = 0
      for (i in 1:n){
        somme = somme + Y[i]*base_trigonometrique(X[i], j)
      }
      return(somme/n)
    }
    
    # MATRICE DE PROJECTION G 
    
    G = function(j, k){
      somme = 0
      for(i in 1:n){
        somme = somme + base_trigonometrique(X[i], j)*base_trigonometrique(X[i], k)
      }
      return(somme/n)
    }
    
    # L'ESTIMATEUR DE TETA
    
    Teta_chapeau = function(j, N){
      TETA = solve(G(j, j), Y_base_trigonometrique(j))
      return(TETA)
    }
    
    # L'ESTIMATEUR DE f
    
    f_chapeau = function(X, N){
      somme = 0
      for (i in 1:N){
        somme = somme + Teta_chapeau(i, N)*base_trigonometrique(X, i)
      }
      return(somme)
    }
    
    # REPRESENTATION GRAPHIQUE DE L'ESTIMATEUR
    
    t = seq(from = 0, to = 1, length = 1000)
    
    data = data.frame(t = t, f = f(t), f_chapeau = sapply(t, f_chapeau, N = N)) 
    
    ggplot(data, aes(x = t)) +
      
      geom_line(aes(y = f), color = "blue", linetype = "solid") +
      
      geom_line(aes(y = f_chapeau), color = "red", linetype = "dashed") +
      
      labs(title = "Fonction f_chapeau", x = "t", y = "Valeur") +
      
      theme_minimal()
    
