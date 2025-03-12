library(shiny)

# Function to extract HOMO and LUMO energies from Gaussian log file
extract_homo_lumo <- function(logfile) {
  lines <- readLines(logfile)
  homo_values <- c()
  lumo_values <- c()
  
  for (i in seq_along(lines)) {
    if (grepl("Alpha  occ. eigenvalues", lines[i])) {
      values <- as.numeric(na.omit(as.numeric(unlist(strsplit(trimws(lines[i]), "[[:space:]]+")))))
      homo_values <- c(homo_values, values)
    }
    if (grepl("Alpha virt. eigenvalues", lines[i])) {
      values <- as.numeric(na.omit(as.numeric(unlist(strsplit(trimws(lines[i]), "[[:space:]]+")))))
      lumo_values <- c(lumo_values, values)
    }
  }
  
  if (length(homo_values) > 0) {
    homo <- tail(homo_values, 1)
  } else {
    homo <- NA
  }
  if (length(lumo_values) > 0) {
    lumo <- head(lumo_values, 1)
  } else {
    lumo <- NA
  }
  
  return(list(HOMO = homo, LUMO = lumo))
}

# Function to extract vibrational frequencies and thermodynamic properties
extract_vibrations <- function(logfile) {
  lines <- readLines(logfile)
  freqs <- c()
  enthalpy <- NA
  gibbs <- NA
  
  for (i in seq_along(lines)) {
    if (grepl("Frequencies --", lines[i])) {
      values <- as.numeric(unlist(strsplit(trimws(lines[i]), "[[:space:]]+"))[3:5])
      freqs <- c(freqs, values)
    }
    if (grepl("Enthalpy=", lines[i])) {
      enthalpy <- as.numeric(unlist(strsplit(trimws(lines[i]), "[[:space:]]+"))[2])
    }
    if (grepl("Gibbs Free Energy=", lines[i])) {
      gibbs <- as.numeric(unlist(strsplit(trimws(lines[i]), "[[:space:]]+"))[4])
    }
  }
  
  return(list(Frequencies = freqs, Enthalpy = enthalpy, Gibbs = gibbs))
}

# Function to compute derived properties
calculate_properties <- function(homo, lumo) {
  Eg <- lumo - homo  # Energy gap
  IE <- -homo        # Ionization potential
  EA <- -lumo        # Electron affinity
  chi <- (IE + EA) / 2  # Electronegativity
  mu <- -chi         # Electrochemical potential
  eta <- Eg / 2      # Hardness
  sigma <- 1 / eta   # Softness
  omega <- chi^2 / (2 * eta)  # Electrophilicity Index
  
  return(list(Eg = Eg, IE = IE, EA = EA, chi = chi, mu = mu, eta = eta, sigma = sigma, omega = omega))
}

# Shiny UI
ui <- fluidPage(
  titlePanel("Gaussian Log File Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload Gaussian Log File", accept = ".log")
    ),
    mainPanel(
      tableOutput("results"),
      tableOutput("vibrations")
    )
  )
)

# Shiny Server
server <- function(input, output) {
  output$results <- renderTable({
    req(input$file)
    
    log_data <- extract_homo_lumo(input$file$datapath)
    
    if (is.na(log_data$HOMO) || is.na(log_data$LUMO)) {
      return(data.frame(Error = "HOMO/LUMO values not found"))
    }
    
    properties <- calculate_properties(log_data$HOMO, log_data$LUMO)
    
    results <- data.frame(
      Property = c("HOMO (eV)", "LUMO (eV)", "Energy Gap (eV)", "Ionization Potential (IE) (eV)", 
                   "Electron Affinity (EA) (eV)", "Electronegativity (𝜒) (eV)", 
                   "Electrochemical Potential (𝜇) (eV)", "Hardness (𝜂) (eV)", "Softness (𝜎) (eV)", "Electrophilicity (ω) (eV)"),
      Value = c(log_data$HOMO, log_data$LUMO, properties$Eg, properties$IE, properties$EA, 
                properties$chi, properties$mu, properties$eta, properties$sigma, properties$omega)
    )
    return(results)
  })
  
  output$vibrations <- renderTable({
    req(input$file)
    vib_data <- extract_vibrations(input$file$datapath)
    
    results <- data.frame(
      Property = c("First 3 Vibrational Frequencies (cm-1)", "Enthalpy (Hartree)", "Gibbs Free Energy (Hartree)"),
      Value = c(paste(vib_data$Frequencies, collapse = ", "), vib_data$Enthalpy, vib_data$Gibbs)
    )
    return(results)
  })
}

shinyApp(ui, server)
