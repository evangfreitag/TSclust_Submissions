# Updated diss.CID function

diss.CID <- function (x, y) {
    #.ts.sanity.check(x, y)
    CED.x = sqrt(sum(diff(x)^2))
    CED.y = sqrt(sum(diff(y)^2))
    denom <- min(CED.x, CED.y)
    if(denom == 0) {stop("Cannot divide by zero: A series exists that has complexity zero.")}
    else if(denom > 0) {
    CF = max(CED.x, CED.y)/denom
    CF * dist(rbind(x, y))
    }
}
