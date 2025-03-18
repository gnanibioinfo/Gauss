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
  
  homo <- ifelse(length(homo_values) > 0, tail(homo_values, 1), NA)
  lumo <- ifelse(length(lumo_values) > 0, head(lumo_values, 1), NA)
  
  return(list(HOMO = homo, LUMO = lumo))
}

# Function to extract dipole moment and polarizability
dipole_polarizability <- function(logfile) {
  lines <- readLines(logfile)
  dipole <- NA
  polarizability <- NA
  
  for (i in seq_along(lines)) {
    if (grepl("Dipole moment", lines[i])) {
      dipole <- as.numeric(unlist(strsplit(trimws(lines[i]), "[[:space:]]+"))[4])
    }
    if (grepl("Isotropic polarizability", lines[i])) {
      polarizability <- as.numeric(unlist(strsplit(trimws(lines[i]), "[[:space:]]+"))[4])
    }
  }
  
  return(list(Dipole = dipole, Polarizability = polarizability))
}

# Function to extract thermodynamic properties
extract_thermo <- function(logfile) {
  lines <- readLines(logfile)
  enthalpy <- NA
  gibbs <- NA
  entropy <- NA
  
  for (i in seq_along(lines)) {
    if (grepl("Enthalpy=", lines[i])) {
      enthalpy <- as.numeric(unlist(strsplit(trimws(lines[i]), "[[:space:]]+"))[2])
    }
    if (grepl("Gibbs Free Energy=", lines[i])) {
      gibbs <- as.numeric(unlist(strsplit(trimws(lines[i]), "[[:space:]]+"))[4])
    }
    if (grepl("Entropy=", lines[i])) {
      entropy <- as.numeric(unlist(strsplit(trimws(lines[i]), "[[:space:]]+"))[2])
    }
  }
  
  return(list(Enthalpy = enthalpy, Gibbs = gibbs, Entropy = entropy))
}

# Function to compute derived properties
calculate_properties <- function(homo, lumo) {
  Eg <- lumo - homo
  IE <- -homo
  EA <- -lumo
  chi <- (IE + EA) / 2
  mu <- -chi
  eta <- Eg / 2
  sigma <- 1 / eta
  omega <- chi^2 / (2 * eta)
  
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
      tableOutput("thermo"),
      tableOutput("dipole")
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
  
  output$thermo <- renderTable({
    req(input$file)
    thermo_data <- extract_thermo(input$file$datapath)
    
    results <- data.frame(
      Property = c("Enthalpy (Hartree)", "Gibbs Free Energy (Hartree)", "Entropy (cal/mol-K)"),
      Value = c(thermo_data$Enthalpy, thermo_data$Gibbs, thermo_data$Entropy)
    )
    return(results)
  })
  
  output$dipole <- renderTable({
    req(input$file)
    dipole_data <- dipole_polarizability(input$file$datapath)
    
    results <- data.frame(
      Property = c("Dipole Moment (Debye)", "Polarizability (a.u.)"),
      Value = c(dipole_data$Dipole, dipole_data$Polarizability)
    )
    return(results)
  })
}

shinyApp(ui, server)
