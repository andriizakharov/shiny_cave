### Presets for input values
### These three are from the paper, but we need the ones from other papers

preset1 <- c(
    0, -14, 100, 200, 25, 0, # variable X
    0, -5, 100, 200, 25, 0, # variable Y
    -3, 0, 0, 0.5, # inter-variable
    10, 20 # timepoints and simulation sample
)

preset2 <- c(
    0, -14, 100, 200, 25, 0, # variable X
    0, -2.5, 100, 200, 25, 0, # variable Y
    -3, 0, 0, 0.5, # inter-variable
    10, 20 # timepoints and simulation sample
)

preset3 <- c(
    0, -14, 100, 200, 25, 0, # variable X
    0, -5, 100, 200, 25, 0, # variable Y
    60, 0, 0, 0.5, # inter-variable
    10, 20 # timepoints and simulation sample
)

# TBD?
# names(preset1) <- c("icept", "mean", ... )
# names(preset2) <- names(preset1)
# names(preset3) <- names(preset1)