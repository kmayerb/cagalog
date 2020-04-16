require(ggplot2)
require(dplyr)
require(tidyr)
require(ggrepel)

volcanic <- function(df, var){

    # SUBSET TO A SPECIFIC VAR
    df2 <- df %>% dplyr::filter(parameter == var )%>% 
        tidyr::spread(key = type , value = value) %>% 
        dplyr::mutate(log10pv = -1*log10(p_value)) %>% 
        dplyr::arrange(desc(log10pv)) 


    nrows = dim(df2)[1]
    topn = 100
    set.seed(1)
    random_sampling = dplyr::sample_frac(df2, size= .1)
    top_hits = df2[1:topn,]
    df2pres = rbind(top_hits, random_sampling)

    df2pres[['color']] <- rep("sample", dim(df2pres)[1])
    df2pres[['color']][1:topn] <- "top hits"

    gg <- ggplot(df2pres, aes(x = estimate, y = log10pv, col = color, size = color)) + geom_point(alpha = .1) + 
    scale_color_manual(values = c("black", "red")) +
    scale_size_manual(values = c(.1,2)) +
    theme_classic() + 
    theme(legend.position = "none") + 
    geom_text_repel(data = top_hits[0:10,], aes(label = CAG, color = NULL), size = 4 ) + 
    ggtitle(paste0(var , ": Top ", topn ,"\nHits and Random 10% of Data")) + 
    ylab("-1*log10(P Value)") 
            
    return(gg)
}





fn = '/Volumes/LaCie/Users/kmayerbl/gscf/stats/corncob.results.csv' 
df = readr::read_csv(fn)
variables = unique(df$parameter)
var = "mu.cf_statusControl"
x = volcanic(df = df , var = "mu.cf_statusControl")

ggs = purrr::map(variables[0:12], ~volcanic(df = df , var = .x))

i = 0
for (v in variables[0:12]){
    print(v)
    i = i + 1
    #v = stringr::str_replace(v, pattern= ".","_")
    pdf(paste0("volcano_figs/", v,".volc.pdf"), height = 8, width = 8)
    print(ggs[i])
    dev.off()
}

i = 0
pdf(paste0("volcano_figs/all_volcanos.pdf"), height = 8, width = 8)
for (v in variables[0:12]){
    print(v)
    i = i + 1
    #v = stringr::str_replace(v, pattern= ".","_")
    print(ggs[i])
}
dev.off()
