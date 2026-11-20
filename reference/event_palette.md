# Default fill colours for splicelogic event types

One pastel colour per `event_type` returned by the `find_*()` functions,
chosen together with the soft blue and red
[`plot_transcripts()`](https://thelovelab.github.io/splicelogic/reference/plot_transcripts.md)
uses for the up- and down-regulated transcripts, so the thirteen stay as
far apart as pastels allow, including under protanopia and deuteranopia.
Pastels leave less room between colours than saturated ones, which is
why event exons are also outlined and labelled with their type.

## Usage

``` r
event_palette()
```

## Value

A character vector of colours, named by event type.

## Examples

``` r
event_palette()
#>        se        ie       cse       cie       mxe        ri        si      a5ss 
#> "#FFB38A" "#74B9B7" "#E87C9C" "#B784C4" "#AAA7F1" "#DDE699" "#62F9EA" "#7BE8FE" 
#>      a3ss      atss      ates 
#> "#FFC4D6" "#99D6A8" "#81bc5c" 
```
