

# Scratch work---nothing special here I don't think.

benin.dta = subset( benin.dta, Treatment != "Clientelist" )
benin.dta

tapply( benin.dta$VoteShare, benin.dta$Treatment, mean )

rks = rank( benin.dta$VoteShare )
rks
plot( rks ~ I( as.numeric( benin.dta$Treatment ) / 2), xlim=c(0,6) )


sr = StephensonRanks( -benin.dta$VoteShare, 6 )
sr

stripchart( sr ~ benin.dta$Treatment, method="stack", pch=19,offset=2/3 )



SR = function( y, tx , size ) {
    rk = StephensonRanks( y, size )
    mns = tapply( rk, tx, sum )
   # stripchart( rk ~ tx, method="stack", pch=19,offset=2/3 )

    mns[ length(mns) ] #- mns[2]
}


do.test = function( y, tx, size=6 ) {

    tobs = SR( y, tx, size )
    ptest = replicate( 1000, {
        tstar = sample( tx )
        SR( y, tstar, size )
    })
    hist( ptest )
    abline( v=tobs, col="red" )
    mean( ptest > tobs )

}


SR( benin.dta$VoteShare, benin.dta$Treatment, 6 )
SR( -benin.dta$VoteShare, benin.dta$Treatment, 6 )

do.test( benin.dta$VoteShare, benin.dta$Treatment, 4 )
do.test( benin.dta$VoteShare, benin.dta$Treatment, 5 )
do.test( benin.dta$VoteShare, benin.dta$Treatment, 6 )
do.test( benin.dta$VoteShare, benin.dta$Treatment, 7 )
do.test( benin.dta$VoteShare, benin.dta$Treatment, 8 )

do.test( -benin.dta$VoteShare, benin.dta$Treatment, 4 )
do.test( -benin.dta$VoteShare, benin.dta$Treatment, 5 )
do.test( -benin.dta$VoteShare, benin.dta$Treatment, 6 )
do.test( -benin.dta$VoteShare, benin.dta$Treatment, 7 )
do.test( -benin.dta$VoteShare, benin.dta$Treatment, 8 )




y = benin.dta$VoteShare
y[ benin.dta$Treatment == "Policy" ] = y[ benin.dta$Treatment == "Policy" ]  + .18

tobs = SR( y, benin.dta$Treatment, 6 )
tobs
tobs.n = SR( -y, benin.dta$Treatment, 6)
tobs.n

do.test( y, benin.dta$Treatment, 6 )
do.test( -y, benin.dta$Treatment, 6 )


mean( ptest > tobs )




Y1 =  rexp( 40 )
Y2 = rnorm( 40, 1 )
Y2 = Y2 - mean(Y2) + mean(Y1)

t.test( Y1, Y2 )

y = c( Y1, Y2 )
tx = rep( c(0,1), each=20 )
tx
length( y )
length( tx)
SR( y, tx, 6 )

do.test( y, tx, 6 )
do.test(-y, tx, 6 )

boxplot( y ~ tx )
