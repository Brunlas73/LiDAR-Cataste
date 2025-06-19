library(shiny)
library(lidR)
library(rgl)
library(tools)
library(Rvcg)
library(data.table)
library(rlas)
library(dbscan)
library(geometry)
library(ggplot2)
library(dplyr)
library(raster)
library(rgl)
library(data.table)
library(Rvcg)
library(sf)
library(tcltk)
library(terra)
library(tools)
suppressMessages(library(doParallel))
suppressMessages(library(parallel))

# Désactive les fenêtres OpenGL externes de rgl
options(rgl.useNULL = TRUE)

# Permet de charger des fichiers plus gros
options(shiny.maxRequestSize = 300 * 1024^2)
verbatimTextOutput("calc_result")  # à placer dans la partie mainPanel
plotOutput("axes_plot", height = "500px")

ui <- fluidPage(
  titlePanel("Visualisation et Estimation 3D LAS / PLY"),
  sidebarLayout(
    sidebarPanel(
      fileInput("input_file", "Charger un fichier .las ou .ply", accept = c(".las", ".ply")),
      actionButton("go", "Visualisation 3D"),
      br(),
      actionButton("compute", "Calcul")
      # br(),
      # actionButton("go2", "Visualisation plot")
    ),
    
    mainPanel(
      rglwidgetOutput("plot3d_output"),
#      rglwidgetOutput("axes_plot3d"),
      verbatimTextOutput("calc_result")
    )
  )
)

server <- function(input, output, session) {
  hull_points <- reactiveVal(NULL)  
  las_clu <- reactiveVal(NULL)  
  # Réactif principal : lit un LAS ou PLY et retourne un objet LAS
  las_data <- reactive({
    req(input$input_file)
    
    ext <- tolower(file_ext(input$input_file$name))
    path <- input$input_file$datapath
    
    if (ext == "las") {
      las <- readLAS(path)
      validate(need(!is.empty(las), "Fichier LAS vide ou invalide."))
      coords <- as.matrix(las@data[, c("X", "Y", "Z")])
      coords_centered <- scale(coords, center = TRUE, scale = FALSE)  # Centrage
      pca <- prcomp(coords_centered)
      rot_mat <- pca$rotation  # Matrice de rotation 3x3
      rotated_coords <- coords_centered %*% rot_mat
      
      las@data$X <- round(rotated_coords[,1], 3)
      las@data$Y <- round(rotated_coords[,2], 3)
      las@data$Z <- round(rotated_coords[,3], 3)
      return(las)
      
    } else if (ext == "ply") {
      ply <- try(vcgPlyRead(path, updateNormals = FALSE), silent = TRUE)
      if (inherits(ply, "try-error") || is.null(ply$vb)) {
        showNotification("Erreur : fichier PLY invalide.", type = "error")
        return(NULL)
      }
      
      coords <- t(ply$vb[1:3, , drop = FALSE])
      rm(ply)
      gc()
      colnames(coords) <- c("X", "Y", "Z")
      coords <- coords[, c("X", "Z", "Y")]
      colnames(coords) <- c("X", "Y", "Z")
      # Centrage facultatif (mais conseillé si les coordonnées sont grandes)
      coords_centered <- scale(coords, scale = FALSE)
      points_df <- as.data.table(coords_centered)
      setnames(points_df, c("X", "Y", "Z"))  # s'assurer que les noms sont corrects
      #points_df[, Classification := 2L]
      las <- LAS(points_df)
      coords <- as.matrix(las@data[, c("X", "Y", "Z")])
      coords_centered <- scale(coords, center = TRUE, scale = FALSE)  # Centrage
      pca <- prcomp(coords_centered)
      rot_mat <- pca$rotation  # Matrice de rotation 3x3
      rotated_coords <- coords_centered %*% rot_mat
      
      las@data$X <- round(rotated_coords[,1], 3)
      las@data$Y <- round(rotated_coords[,2], 3)
      las@data$Z <- round(rotated_coords[,3], 3)
      return(las)
    }
    
    # Centrage
    # xy <- scale(las[, 1:2], scale = FALSE)  # X, Y seulement
    # # PCA pour trouver les axes propres
      # pca <- prcomp(xy)
    # rotated <- xy %*% pca$rotation  # points tournés selon l'enveloppe
    

    
    showNotification("Extension non supportée.", type = "error")
    return(NULL)
  })
  
  # VISUALISATION 3D
  output$plot3d_output <- renderRglwidget({
    req(input$go)
    las <- las_data()
    req(!is.null(las))
    
    x <- las@data$X
    y <- las@data$Y
    z <- las@data$Z
    validate(need(length(x) > 0, "Aucune coordonnée disponible."))
    
    clear3d()
    lims <- range(c(x, y, z))
    
    plot3d(
      x, y, z,
      col = "darkgreen",
      size = 1,
      xlab = "X", ylab = "Y", zlab = "Z",
      xlim = lims,
      ylim = lims,
      zlim = lims,
      aspect = c(1, 1, 1)
    )
    
    rglwidget()
  })

  output$axes_plot3d <- renderRglwidget({
    req(input$go2)
      las <- las_clu()
      req(!is.null(las))
      
      pts <- las@data[, c("X", "Y", "Z")]

      hull_pts <- hull_points()
      
      if (!is.null(hull)) {
        

      # PCA sur X, Y
      xy_centered <- scale(coords[, 1:2], scale = FALSE)
      pca <- prcomp(xy_centered)
      eigvec <- pca$rotation
      center <- colMeans(hull_pts[, 1:3])
      
      clear3d()
      bg3d(color = "white")
      points3d(hull_pts, color = "black", size = 4)
      points3d(pts, color = "red", size = 1)
      
      # Axes d'origine (gris)
      lines3d(rbind(center, center + c(50, 0, 0)), color = "grey", lwd = 3)  # X
      lines3d(rbind(center, center + c(0, 50, 0)), color = "grey", lwd = 3)  # Y
      lines3d(rbind(center, center + c(0, 0, 50)), color = "grey", lwd = 3)  # Z
      
      # Axes PCA (rouge et bleu) dans le plan XY
      lines3d(rbind(center,
                    center + c(50 * eigvec[1, 1], 50 * eigvec[2, 1], 0)),
              color = "red", lwd = 3)
      lines3d(rbind(center,
                    center + c(50 * eigvec[1, 2], 50 * eigvec[2, 2], 0)),
              color = "blue", lwd = 3)
      }
      rglwidget()
    })
    
  # CALCUL
  calc_result <- eventReactive(input$compute, {
    las <- las_data()
    req(!is.null(las))
    
    # z_vals <- las@data$Z
    # if (length(z_vals) == 0) return("Aucune altitude (Z) disponible.")
    # 
    # stats <- summary(z_vals)
    # paste("Résumé des hauteurs (Z) :", capture.output(stats), collapse = "\n")

        las@data$Classification <- rep(1L, nrow(las@data))
    #lidR::projection(las) <- "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"
    
    dimtu <- 0.4
    
    tiles_sd <- grid_metrics(las, sd(Z), res = dimtu)
    
    # 3. Convertir en data frame avec les coordonnées X et Y
    tiles_df <- as.data.frame(tiles_sd, xy = TRUE)
    
    # 4. Vérifie et renomme les colonnes si nécessaire
    if (ncol(tiles_df) == 3) {
      names(tiles_df) <- c("x", "y", "sd")
    } else {
      stop("Erreur : le raster des tuiles n'a pas été correctement transformé.")
    }
    
    # 5. Filtrer les tuiles plates (seuil modifiable)
    flat_tiles <- subset(tiles_df, sd < 0.05)
    
    # 6. Visualiser les zones plates (optionnel)
    # plot(tiles_df$x, tiles_df$y, col = ifelse(tiles_df$sd < 0.05, "blue", "gray"), pch = 20,
    #      main = "Zones plates détectées", xlab = "X", ylab = "Y")
    
    # 7. Extraire les points dans les tuiles plates (rayon de 0.5 autour de chaque centre de tuile)
    points_sol <- list()
    
    for (i in 1:nrow(flat_tiles)) {
      x <- flat_tiles$x[i]
      y <- flat_tiles$y[i]
      
      # Sélection des points dans la tuile de 1m centrée sur (x, y)
      sub <- filter_poi(las, X >= x - (dimtu/2) & X < x + (dimtu/2) & Y >= y - (dimtu/2) & Y < y + (dimtu/2))
      
      if (!is.null(sub)) {
        # Sélection du point le plus bas (Z minimal)
        zmin <- filter_poi(sub, Z == min(Z, na.rm = TRUE))
        points_sol[[length(points_sol) + 1]] <- zmin
      }
    }
    
    # 8. Fusionner tous les points du sol extraits
    sol_las <- do.call(rbind, points_sol)
    
    # 9. Vérification
    if (is.null(sol_las)) stop("Aucun point du sol n'a pu être extrait.")
    
    # Ajoute une classification "ground" aux points retenus
    sol_las$Classification <- 2L
    
    # 3. Filtrer les points ground trop hauts (ex. au-dessus du 95e percentile)
    z_vals <- sol_las@data$Z
    z_thresh <- quantile(z_vals, 0.90)  # ou 0.90 selon les cas
    
    # Conserver uniquement les points ground en dessous du seuil
    ground_las_clean <- filter_poi(sol_las, Z < z_thresh)
    
    # Définir les bornes complètes avec une petite marge
    xrange <- range(las@data$X)
    yrange <- range(las@data$Y)
    
    # Étendre un peu les bornes
    buffer <- 0.1
    xrange <- xrange + c(-buffer, buffer)
    yrange <- yrange + c(-buffer, buffer)
    
    # Crée un DTM avec ces points
    dtm <- rasterize_terrain(ground_las_clean, algorithm = tin(), res = 0.1, bounds = c(xrange[1], xrange[2], yrange[1], yrange[2]))
    
    # Vérifie visuellement le DTM
    #plot(dtm, main = "DTM basé sur les points plats du sol")
    
    
    # 1. Extraire la couche raster
    dtm_2d <- raster::raster(dtm)
    
    # 2. Convertir en matrice
    z <- raster::as.matrix(dtm_2d)
    
    # 3. Créer les grilles X et Y
    x <- seq(from = xmin(dtm_2d), to = xmax(dtm_2d), length.out = ncol(z))
    y <- seq(from = ymin(dtm_2d), to = ymax(dtm_2d), length.out = nrow(z))
    
    # 4. Créer meshgrid avec expand.grid, puis retransformer en matrice
    X <- matrix(rep(x, each = length(y)), nrow = length(y), byrow = FALSE)
    Y <- matrix(rep(y, times = length(x)), nrow = length(y), byrow = FALSE)
    
    # 5. Afficher le DTM en 3D
    # open3d()
    # surface3d(X, Y, z, color = "tan", back = "lines")
    # points3d(las@data$X, las@data$Y, las@data$Z, col = rgb(1, 0, 1, 0.2), size = 2)
    # points3d(sol_las@data$X, sol_las@data$Y, sol_las@data$Z, col = rgb(0, 0, 1, 0.2), size = 2)
    
    # 8. Visualisation
    #plot(dtm, main = "DTM basé sur zones plates")
    
    # Normaliser le nuage avec ce DTM
    #las_normalized <- normalize_height(las, dtm)
    
    # Exemple : définir un CRS UTM arbitraire, ici EPSG:32633 (UTM zone 33N)
    crs_string <- "+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"
    
    # Appliquer le CRS au LAS
    #projection(las) <- crs_string
    lidR::projection(las) <- crs_string
    #"+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"
    
    # Appliquer le même CRS au DTM si nécessaire
    crs(dtm) <- crs_string
    
    #las <- las_normalized
    cooxyz <- las@data[, c("X", "Y", "Z")]
    z_ground <- terra::extract(dtm, cooxyz[, c("X", "Y")])[,2]
    buffer_height <- 0.1  # 20 cm
    keep <- !is.na(z_ground) & abs(cooxyz[, "Z"] - z_ground) <= buffer_height
    las@data$Classification[keep] <- 2L
    
    ground <- filter_poi(las, Classification == 2)
    non_ground <- filter_poi(las, Classification != 2)
    
    plot_ground_vs_nonground <- function(las) {
      if (!"Classification" %in% colnames(las@data)) {
        stop("⚠️ Le nuage de points n'est pas classifié. Appliquer classify_ground() d'abord.")
      }
      
      # Séparer
      ground <- filter_poi(las, Classification == 2)
      nonground <- filter_poi(las, Classification != 2)
      
      # Ajouter une valeur numérique pour la couleur (0 = sol, 1 = non-sol)
      ground@data$color_code <- 0
      nonground@data$color_code <- 1
      
      # Fusionner
      las_combined <- rbind(ground, nonground)
      
      # Afficher
      # plot(las_combined, color = "color_code", size = 3, bg = "white")
    }# Visualisation sol vs non-sol
    # plot_ground_vs_nonground(las)
    
    # -------- 3. Heuristique automatique du meilleur eps --------
    auto_eps <- function(coords, sample_size = 5000, k = 4, plot = FALSE) {
      library(dbscan)
      
      # Échantillonnage des points pour réduire la charge mémoire
      if (nrow(coords) > sample_size) {
        set.seed(42)
        coords_sample <- coords[sample(1:nrow(coords), sample_size), ]
      } else {
        coords_sample <- coords
      }
      
      # Calcul des distances aux k plus proches voisins
      kNNdist <- kNNdist(coords_sample, k = k)
      kNNdist_sorted <- sort(kNNdist)
      
      # Affichage du graphe pour identifier le "coude"
      if (plot) {
        plot(kNNdist_sorted, type = "l", main = "k-NN Distance Plot", xlab = "Points triés", ylab = paste0(k, "-ème plus proche voisin"))
        grid()
      }
      
      # Estimation heuristique de eps (quantile du coude)
      eps_guess <- quantile(kNNdist_sorted, 0.95)
      return(as.numeric(eps_guess))
    }
    
    
    coords <- as.matrix(non_ground@data[, c("X", "Y", "Z")])
    
    best_eps <- auto_eps(coords)
#    cat("✨ eps optimal détecté :", best_eps, "\n")
    
    set.seed(1)
    sample_idx <- sample(nrow(coords), min(nrow(coords), 100000))
    coords_sample <- coords[sample_idx, ]
    # -------- 4. DBSCAN avec eps optimal --------
    # db <- dbscan(coords_sample, eps = best_eps, minPts = 50, search = "kdtree")
    db <- dbscan(coords, eps = best_eps, minPts = 50, search = "kdtree")
    
    # coords_sample_df <- as.data.frame(coords_sample)
    coords_sample_df <- as.data.frame(coords)
    coords_sample_df$cluster <- db$cluster
    
    # Exclurele b ruit (cluster 0)
    valid_clusters <- coords_sample_df$cluster[coords_sample_df$cluster != 0]
    if (length(valid_clusters) == 0) stop("Aucun cluster détecté")
    
    # Trouver le cluster le plus grand
    cluster_sizes <- table(valid_clusters)
    biggest_cluster_id <- as.integer(names(which.max(cluster_sizes)))
    
    # Extraire les points du plus gros cluster
    biggest_cluster <- coords_sample_df[coords_sample_df$cluster == biggest_cluster_id, ]
    
    # Cluster principal
    X_clust <- biggest_cluster[, 1]
    Y_clust <- biggest_cluster[, 2]
    Z_clust <- biggest_cluster[, 3]
    
    X_las <- cooxyz[, 1]
    Y_las <- cooxyz[, 2]
    Z_las <- cooxyz[, 3]
    
    # Limites identiques pour tous les axes
    lims <- range(c(X_clust, Y_clust, Z_clust))
    
    # # Nuage complet en gris
    # plot3d(X_clust, Y_clust, Z_clust,
    #        col = "blue",
    #        size = 0.1,
    #        type = "s",
    #        xlab = "X", ylab = "Y", zlab = "Z",
    #        xlim = lims,
    #        ylim = lims,
    #        zlim = lims)
    # 
    # points3d(xyz,
    #          col = "red",
    #          size = 0.1)
    # 
    # # Même échelle sur les axes
    # aspect3d(1, 1, 1)
    # 
    las_cluster <- LAS(biggest_cluster)
    summary(las_cluster)
    
    las_clu(las_cluster)
    
    # 1. Extraire les coordonnées
    cooxyz <- las_cluster@data[, c("X", "Y", "Z")]
    
    # 2. Calcul de l'enveloppe convexe (triangles)
    hull <- convhulln(cooxyz, options = "FA")  # "FA" renvoie aussi l'aire et le volume
    # 3. Créer le mesh3d (attention à l’ordre des dimensions)
    vertices <- t(cooxyz)
    faces <- t(hull$hull)  # transposer pour avoir la bonne orientation
    
    mesh <- tmesh3d(
      vertices = vertices,
      indices = faces,
      homogeneous = FALSE
    )
    
    # Extraire les sommets du hull (points uniques)
    hull_pts <- unique(as.data.frame(cooxyz[unique(as.vector(hull$hull)), ]))
    hull_points(hull_pts)
    
        # Centrage
    xy <- scale(hull_pts[, 1:2], scale = FALSE)  # X, Y seulement
    

    # Étendue sur les axes propres
    longueur <- diff(range(xy[, 1]))
    largeur  <- diff(range(xy[, 2]))
    hauteur  <- diff(range(hull_pts[, "Z"]))
        
    # 4. Affichage
    
    # shade3d(mesh, color = "tan")
    
    # 5. Volume (en m³) directement depuis convhulln
    volenv <- round(hull$vol, 2)
#    cat("Volume estimé (enveloppe convexe) :", round(hull$vol, 2), "m³\n")
    
    ######################################
    
    # 1. Génère le DSM du tas de bois à partir du cluster
    dsm <- rasterize_canopy(las_cluster, res = 0.1, algorithm = p2r())
    
    # 2. Recadre le DTM global pour qu'il corresponde à l'emprise du cluster
    dtm_crop <- crop(dtm, dsm)
    dtm_resample <- resample(dtm_crop, dsm, method = "bilinear")
    
    # 3. Calcule la hauteur = DSM - DTM global
    height_raster <- dsm - dtm_resample
    
    # 4. Intègre la hauteur pour obtenir le volume
    volume_m3 <- global(height_raster, fun = "sum", na.rm = TRUE)[[1]] * 0.01  # 0.1m * 0.1m = 0.01 m² par pixel
    
    # 5. Affiche le volume
    #cat("Volume estimé (point par point) :", round(volume_m3, 2), "m³\n")
    volpoi <- round(volume_m3, 2)
    
    #################################################
    
  
  # On capture les lignes comme si on faisait cat()
  lignes <- c(
    paste("Volume estimé (enveloppe convexe) :", volenv),
    paste("Largeur  :", round(largeur, 2), "m"),
    paste("Longueur :", round(longueur, 2), "m"),
    paste("Hauteur  :", round(hauteur, 2), "m"),
    paste("Volume estimé (point par point) :", volpoi)
  )
  
  paste(lignes, collapse = "\n")
})

  output$calc_result <- renderText({
    calc_result()
  })
}

shinyApp(ui, server)
