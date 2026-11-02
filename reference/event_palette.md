# Default fill colours for splicelogic event types

One pastel colour per `event_type` returned by the `find_*()` functions,
chosen together with the soft blue and red
[`plot_transcripts()`](https://thelovelab.github.io/splicelogic/reference/plot_transcripts.md)
uses for the up- and down-regulated transcripts, so the eleven stay as
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
#>        se        ie       mxe        ri        si      a5ss      a3ss      atss 
#> "#FFB38A" "#74B9B7" "#AAA7F1" "#DDE699" "#62F9EA" "#7BE8FE" "#FFC4D6" "#99D6A8" 
#>      ates 
#> "#81bc5c" 
```
