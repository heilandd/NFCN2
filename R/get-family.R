#' @title getPoints
#' @author Dieter Henrik Heiland
#' @description getPoints
#' @inherit 
#' @param coordinates A data.frame will rownames as identifier and coordinates names as x and y
#' @return 
#' @examples 
#' 
#' @export

getPoints <- function(object, reduction="umap"){
  
  if(class(object)=="data.frame"){ 
    message(" Input is a data.frame ")
    names(object)[1:2] <- c("x", "y")
    coordinates <- object[,1:2]
    
    }
  if(class(object)=="Seurat"){
    message("Input a Seurat object")
    
    if(base::any(names(object@reductions)==reduction )){
      
      coordinates <- object@reductions[[reduction]]@cell.embeddings[, 1:2] %>% as.data.frame()
      names(coordinates)[1:2] <- c("x", "y")
    }
    
  }
  
  
  mat_out <- coordinates
  
  ui <- shiny::fluidPage(
    shiny::titlePanel("Choose your Cells of Interest for further analysis"),
    
    # Sidebar layout with input and output definitions ----
    shiny::sidebarLayout(
      
      # Sidebar panel for inputs ----
      shiny::sidebarPanel(
        # done button
        shiny::actionButton("choose_toggle", "Choose/unchoose"),
        # clear button
        shiny::actionButton("reset", "Clear"),
        # done button
        shiny::actionButton("done", "Done"),
        shiny::h3("Instructions:"),
        shiny::tags$ol(
          shiny::tags$li("Highlight Cells from your Input DimRed Plot by ",
                         "clicking and dragging."),
          shiny::tags$li("Click the 'Choose/unchoose' button."),
          
          shiny::tags$li("Click 'Done'.")
        ),
        shiny::h4("Details:"),
        shiny::tags$ul(
          shiny::tags$li(paste("Output of the function is the selected",
                               "Cell names")),
          shiny::tags$li("To start over, click 'Clear'"),
          shiny::tags$li(paste("You can also choose/unchoose specific Cells",
                               "by clicking on them directly"))
        )
        
      ),
      
      # Main panel for displaying outputs ----
      shiny::mainPanel(
        shiny::plotOutput("plot1", height = 350,
                          click = "plot1_click",
                          brush = shiny::brushOpts(id = "plot1_brush"))
      )
    )
  )
  
  
  
  
  
  
  
  
  server <- function(input, output, session) {
    
    ica_space_df=mat_out
    vals <- shiny::reactiveValues(
      keeprows = rep(TRUE, nrow(ica_space_df)))
    
    output$plot1 <- shiny::renderPlot({
      keep    <- ica_space_df[ vals$keeprows, , drop = FALSE]
      exclude <- ica_space_df[!vals$keeprows, , drop = FALSE]
      
      ggplot(keep, aes(x, y)) +
        geom_point(data = ica_space_df, aes(x=x, y = y),
                   size = .5, color = "gray", alpha = .3) +
        geom_point(alpha = .7) +
        geom_point(data = exclude, shape = 21, fill = "red", color = "red") +
        labs(x="Component 1", y="Component 2")+
        theme_void()
    }, height = function() {
      session$clientData$output_plot1_width
    })
    
    # Toggle points that are clicked
    shiny::observeEvent(input$plot1_click, {
      res <- shiny::nearPoints(ica_space_df,
                               xvar = "x",
                               yvar = "y",
                               input$plot1_click,
                               allRows = TRUE)
      
      vals$keeprows <- xor(vals$keeprows, res$selected_)
    })
    
    # Toggle points that are brushed, when button is clicked
    shiny::observeEvent(input$choose_toggle, {
      res <- shiny::brushedPoints(ica_space_df, input$plot1_brush,
                                  xvar = "x",
                                  yvar = "y",
                                  allRows = TRUE)
      
      vals$keeprows <- xor(vals$keeprows, res$selected_)
    })
    
    # Reset all points
    shiny::observeEvent(input$reset, {
      vals$keeprows <- rep(TRUE, nrow(ica_space_df))
    })
    
    shiny::observeEvent(input$done, {
      shiny::stopApp(vals$keeprows)
    })
    
    
    
  }
  sel <- shiny::runApp(shiny::shinyApp(ui, server))
  
  return(rownames(mat_out)[!sel])
}


#' @title getTargets
#' @author Dieter Henrik Heiland
#' @description getTargets
#' @inherit 
#' @param coordinates A data.frame will rownames as identifier and coordinates names as x and y
#' @return 
#' @examples 
#' 
#' @export
getTargets <- function(object){object@assays$NFCN$target$cells}

#' @title getSource
#' @author Dieter Henrik Heiland
#' @description getSource
#' @inherit 
#' @param coordinates A data.frame will rownames as identifier and coordinates names as x and y
#' @return 
#' @examples 
#' 
#' @export
getSource <- function(object){object@assays$NFCN$source$cells}


#' @title getCleaned
#' @author Dieter Henrik Heiland
#' @description getCleaned
#' @inherit 
#' @param coordinates A data.frame will rownames as identifier and coordinates names as x and y
#' @return 
#' @examples 
#' 
#' @export

getCleaned <- function(df, feat, q=0.1){
  
  # clean data by crop the quantils
  
  df <- 
    df %>%
    dplyr::filter(!!sym(feat)>stats::quantile(!!sym(feat), q)) %>% 
    dplyr::filter(!!sym(feat)<stats::quantile(!!sym(feat), 1-q))
  
  return(df)
}








