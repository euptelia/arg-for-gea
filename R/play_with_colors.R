# Choosing colors
# https://www.thinkingondata.com/something-about-viridis-library/
# https://cran.r-project.org/web/packages/scales/scales.pdf
#The full name of a viridis palette: "viridis", "magma", "inferno", or "plasma".
library(scales)

# Given number of categories
n_colors = 20
show_col(viridis_pal()(n_colors))


col_numeric("viridis", domain = NULL)(seq(n_colors))
rev(col_numeric("viridis", domain = NULL)(seq(n_colors)))
show_col(col_numeric("viridis", domain = NULL)(seq(n_colors)))

col_numeric("magma", domain = NULL)(seq(n_colors))
rev(col_numeric("magma", domain = NULL)(seq(n_colors)))
show_col(col_numeric("magma", domain = NULL)(seq(n_colors)))

col_numeric("inferno", domain = NULL)(seq(n_colors))
show_col(col_numeric("inferno", domain = NULL)(seq(n_colors)))

col_numeric("plasma", domain = NULL)(seq(n_colors))
show_col(col_numeric("plasma", domain = NULL)(seq(n_colors)))

# Exponential distribution, mapped continuously
show_col(col_numeric("Blues", domain = NULL)(sort(rexp(16))))

# custom data
test = c(0.0,  19997.8,  39995.6,  59993.4,  79991.2,  99989.0,
         119986.8, 139984.6, 159982.4, 179980.2, 199978.0)
col_numeric("magma", domain = NULL)(test)
show_col(col_numeric("magma", domain = NULL)(test))


