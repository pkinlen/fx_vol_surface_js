The javascript fx vol surface, by Philip Kinlen

To see the script in action look at the following google sheet (template)
https://docs.google.com/spreadsheet/ccc?key=0Aq5Pj_il0hRHdHRtc3YtbTJGVGtKWmFkdlI5NkxuRmc&usp=sharing

The javascript shows how an FX vol surface can be constructed from market inputs
( ATMF vols, risk-reversals and strangle premiums )

Suppose you have the implied vol at a given delta, the maths gets a bit interesting
when you realize that you would normally need the vol to work out the delta.

The important thing is to solve in a clever way so that the surface can be efficiently built.

Here is a bit of background:

For a currency pair, say A and B
we can look at the vanilla options market on that FX rate and see the prices of options 
for a range of strike and maturities.
For each option, since we have the price, strike and maturity,
we can use the Black-Scholes formula to back-out the volatility implied by the price
Then if we plot those vols as a function of strike and maturities,
we have our implied vol surface.



