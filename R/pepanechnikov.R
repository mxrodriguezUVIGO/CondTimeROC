pepanechnikov <-
function(x) {
	res <- sapply(x, function(x) {
		if(x < -1)
			res = 0
		else if (x > 1)
			res = 1
		else
			res = (3/4)*(x - x^3/3) + 1/2
	})
	res
}
