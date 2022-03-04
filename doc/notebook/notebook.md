Notebook
================
Roland Knapp
2022-03-03

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
    ## ✓ tibble  3.1.6     ✓ dplyr   1.0.8
    ## ✓ tidyr   1.2.0     ✓ stringr 1.4.0
    ## ✓ readr   2.1.2     ✓ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

## Predictors of frog survival following translocation

### Dataset structure

-   site_id: site to which frogs were translocated
-   elevation
-   release_date
-   release_day
-   release_year
-   siteid_donor
-   order: first translocation to site = 0, subsequent translocation = 1
-   type: translocation conducted on foot versus by helicopter
-   shore: suitability of shore/bank habitat for frog overwintering
-   pit_tag_ref
-   survival: individual-level estimated survival 1 year following
    translocation, based on cmr surveys conducted for at least two years
    post-translocation
-   sex
-   length: at release
-   weight: at release
-   condition: index based on length and weight
-   swab_id
-   bd_load: at release

Maximum depth and surface area would seem useful to include, but 74976
is a stream/meadow habitat and as such depth and area at that site are
likely not comparable to depth and area at lake habitats.

Frog condition might be problematic because size and residuals are
correlated. Specifically, the largest frogs have the largest positive
residuals from a weight:length regression line, even after log10
transformation. As such, length and condition are probably not
independent.

``` r
read_csv(here::here("data", "clean", "frog_translocation.csv")) %>% 
ggplot(aes(x = log10(length), y = log10(weight))) +
  geom_point() +
  geom_smooth(method = lm)
  plot
```

    ## Rows: 779 Columns: 13
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (4): type, pit_tag_ref, sex, swab_id
    ## dbl  (8): site_id, collect_siteid, length, weight, bd_load, median_survival,...
    ## date (1): release_date
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `geom_smooth()` using formula 'y ~ x'

### Model structure

survival \~ elevation + siteid_donor + release_day + release_year +
order + shore + sex + length + condition + log(bdload) \| site_id

Group-level effects: In addition to the inclusion of site_id as a
group-level effect, may need to include release_year, perhaps nested
within site_id? Given that some sites only have a single release_year,
will this cause problems?
