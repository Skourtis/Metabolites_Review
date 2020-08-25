#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
if (!require("pacman")) install.packages("pacman")
pacman::p_load(shiny, tidyverse, ggplot2, shinydashboard,
               plotly, htmlwidgets, heatmaply, RColorBrewer, visNetwork, shinythemes)
load("./Trial_network.RData")
Change_cell_type_names <- function(x){
    case_when(
    x == "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE" ~ "BLOOD",
    x == "LARGE_INTESTINE" ~ "L_INTESTINE",
    x == "SMALL_INTESTINE" ~ "S_INTESTINE",
    x == "UPPER_AERODIGESTIVE_TRACT" ~ "AERODIGESTIVE",
    x == "GASTROINTESTINAL_TRACT" ~ "GASTROINTESTINAL",
    x == "CENTRAL_NERVOUS_SYSTEM" ~ "CNS",
    x == "AUTONOMIC_GANGLIA" ~ "GANGLIA",
    x == "BILIARY_TRACT" ~ "BILIARY",
    TRUE ~ as.character(x)
)}

load("./CCLE_enzymes_Deepmap_pivoted.Rdata") 
load("./HeatMap.Rdata") 
load("./interactive_ggplot_CCLE_enzymes_test.RData")
Origin <- read_tsv("./Cell_lines_annotations_20181226.txt")
CCLE_enzymes_pivoted$Cell_type <- Change_cell_type_names(CCLE_enzymes_pivoted$Cell_type)
deep_map_pivoted$Cell_type <- Change_cell_type_names(deep_map_pivoted$Cell_type)

# Define UI for application that draws a histogram
header <- dashboardHeader(title = "Sdelci")
sidebar <- dashboardSidebar(
    sidebarMenu(
        menuItem("HeatMap", tabName = "Heatmap"),
        menuItem("Enzymes", tabName = "enzymes"),
        menuItem("Metabolites", tabName = "metabolites"),
        menuItem("Folate Pathway", tabName = "pathway"),
        sidebarSearchForm(textId = "searchText", buttonId = "searchMetabolite", label = "KEGG ID...."),
        sidebarSearchForm(textId = "searchText", buttonId = "searchEnzyme", label = "Gene Name...")
    )
)
body <- dashboardBody({
    tabItems(
        tabItem(tabName = "Heatmap",
                h2("HeatMap"),
    fluidPage(
        fluidRow(
            box(width = 12,
                title = "Metabolites",
                collapsible = T,
                color = "green", ribbon = TRUE, title_side = "top right",
                column(width = 12, plotOutput("HeatMap", height = "1200px")))))),
    tabItem(tabName = "enzymes",
            h2("Enzymes"),
        fluidRow(
            selectInput("SelectedGene", "Select Enzyme:", unique(CCLE_enzymes_pivoted$Gene_Symbol)),
        fluidRow( plotlyOutput("CCLE_enzymes")),
        fluidRow( plotlyOutput("Depmap")))),
    tabItem(tabName = "metabolites",
            h2("Metabolites"),
            fluidRow(
                selectInput("SelectedMetabolite", "Select Metabolite:", unique(folate_cell_line_data_test$Metabolite)),
                fluidRow( plotlyOutput("Metabolites"))),
            fluidRow(
                DT::dataTableOutput("dt")
                )),
    tabItem(tabName = "pathway",
            h2("Folate Pathway"),
        fluidRow(visNetworkOutput("network", height = 1000)) 
            ))
})




# Define server logic required to draw a histogram
server <- function(input, output) {

    output$CCLE_enzymes <- renderPlotly({
        ggplotly_enzymes <-CCLE_enzymes_pivoted[CCLE_enzymes_pivoted$Gene_Symbol == input$SelectedGene,] %>% 
            group_by(Cell_type)%>%
            ggplot(aes(x = Cell_type, y = Abundance, color= Cell_type, group = Cell_type, 
                       text =paste(
                           "Cell Line: ", Cell_line, "\n",
                           "Protein Abundance: ", Abundance, "\n",
                           "Origin: ", "", "\n",
                           "Mutations", "", "\n",
                           sep = ""))) +
            geom_boxplot()+ 
            scale_x_discrete(guide = guide_axis(n.dodge=2))+
            theme(axis.text=element_text(size=7), legend.position='none')
        ggplotly(ggplotly_enzymes, tooltip = "text")  
    })
    
    output$Depmap <- renderPlotly({
       ggplotly_Depmap <-  deep_map_pivoted[deep_map_pivoted$Gene_name == input$SelectedGene,] %>% 
           left_join(CCLE_Mut_profile, by = "Cell_line") %>% 
           left_join(Origin[,1:4], by = c("Cell_line" = "CCLE_ID")) %>%
            group_by(Cell_type)%>%
            ggplot(aes(x = Cell_type, y = Sensitivity, color= Cell_type, group = Cell_type,
                       text =paste(
                           "Cell Line: ", Cell_line, "\n",
                           "CRISPRi Sensitivity: ", Sensitivity, "\n",
                           "Origin: ", Pathology, "\n",
                           "Mutations:", Mutations_Cosmic, "\n",
                           sep = ""))) +
           geom_violin()+
           geom_jitter(aes(shape = Pathology), width = 0.2)+ 
            scale_shape_discrete(solid=F)+
            scale_x_discrete(guide = guide_axis(n.dodge=2))+
            theme(axis.text.x = element_text(size = 7),legend.position='none')
            ggplotly(ggplotly_Depmap, tooltip = "text")  
    })
    
    output$Metabolites <- renderPlotly({
        ggplotly_Metabolites <-  folate_cell_line_data_test[folate_cell_line_data_test$Metabolite == input$SelectedMetabolite,] %>% 
            group_by(Organ)%>%
            ggplot(aes(x = Organ, y = Abundance, color= Organ, group = Organ,
                       text =paste(
                           "Cell Line: ", Organ, "\n",
                           "Abundance: ", Abundance, "\n",
                           sep = ""))) +
            geom_violin()+
            geom_jitter(width = 0.2)+ 
            scale_shape_discrete(solid=F)+
            scale_x_discrete(guide = guide_axis(n.dodge=2))+
            labs(title= input$SelectedMetabolite) +
            theme(axis.text.x = element_text(size = 7),legend.position='none', plot.title = element_text(hjust = 0.5)) 
        ggplotly(ggplotly_Metabolites, tooltip = "text")  
    })
    
    output$HeatMap <- renderPlot({
        heatmap
    })
    output$Enzymes <- renderPlotly({
        Protein_heatmap
            })
    
    output$network <- renderVisNetwork({
        visNetwork(nodes,edges_enzymes, height = 1000, width = "100%") %>% 
            visIgraphLayout(smooth = T)%>% 
            visNodes(
                shape = "box",
                shadow = list(enabled = TRUE, size = 10)
            ) %>%
            visEdges(smooth = list(type = "dynamic"), color = list(highlight = "#C62F4B",opacity = 0.35, border = "black"), 
                     font = list(align = "middle")) %>%
            visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T),
                       selectedBy = "label") %>%
            visPhysics(stabilization = FALSE,solver = "forceAtlas2Based", 
                       forceAtlas2Based = list(gravitationalConstant = -90))
     })
   observe({
        if(input$searchMetabolite > 0){
            isolate({
                print(input$searchText)
                current_node <- nodes[grep(input$searchText, nodes$id), "id"]
                print(current_node)
                visNetworkProxy("network") %>% visSelectNodes(id  = current_node)
            })
        }
    })
    
    output$dt <- DT::renderDataTable({
        DT::datatable(subset(Enzymes_Metabolites, 
                             KEGG == Metabolites_mapped$KEGG[Metabolites_mapped$Query == input$SelectedMetabolite])
                      , select="single")})
    
    observe({
        if(input$searchEnzyme > 0){
        isolate({
            # current_edge <- edges_enzymes[grep(input$searchEnzyme, edges_enzymes$label), "label"]
            # print(current_edge)
            visNetworkProxy("network") %>% visSelectEdges(id  = input$searchEnzyme)
        })
        }
    })
}

# Run the application 
shinyApp(ui = dashboardPage(
    dashboardHeader(title = "Sdelci Lab"),
    sidebar,
    body
), server = server)
