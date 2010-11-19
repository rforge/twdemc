plotkde.2d <- function (fhat, display = "slice", cont = c(25, 50, 75), abs.cont, 
	approx.cont = FALSE, xlab, ylab, zlab = "Density function", 
	cex = 1, pch = 1, labcex, add = FALSE, drawpoints = FALSE, 
	drawlabels = TRUE, theta = -30, phi = 40, d = 4, ptcol = "blue", 
	col, lwd = 1, ...) 
{
	disp1 <- substr(display, 1, 1)
	if (!is.list(fhat$eval.points)) 
		stop("Need a grid of density estimates")
	if (missing(xlab)) 
		xlab <- fhat$names[1]
	if (missing(ylab)) 
		ylab <- fhat$names[2]
	if (missing(labcex)) 
		labcex <- 1
	if (missing(approx.cont)) 
		approx.cont <- (nrow(fhat$x) > 2000)
	if (disp1 == "p") 
		plotret <- persp(fhat$eval.points[[1]], fhat$eval.points[[2]], 
			fhat$estimate, theta = theta, phi = phi, d = d, xlab = xlab, 
			ylab = ylab, zlab = zlab, ...)
	else if (disp1 == "s") {
		if (!add) 
			plot(fhat$x[, 1], fhat$x[, 2], type = "n", xlab = xlab, 
				ylab = ylab, ...)
		if (missing(abs.cont)) {
			if (!is.null(fhat$cont)) {
				cont.ind <- rep(FALSE, length(fhat$cont))
				for (j in 1:length(cont)) cont.ind[which(cont[j] == 
								as.numeric(unlist(strsplit(names(fhat$cont), 
											"%"))))] <- TRUE
				if (all(!cont.ind)) 
					hts <- contourLevels(fhat, prob = (100 - cont)/100, 
						approx = approx.cont)
				else hts <- fhat$cont[cont.ind]
			}
			else hts <- contourLevels(fhat, prob = (100 - cont)/100, 
					approx = approx.cont)
		}
		else hts <- abs.cont
		hts <- sort(hts)
		if (missing(col)) 
			col <- 1
		if (length(col) < length(hts)) 
			col <- rep(col, times = length(hts))
		for (i in 1:length(hts)) {
			if (missing(abs.cont)) 
				scale <- cont[i]/hts[i]
			else scale <- 1
			if (hts[i] > 0) 
				contour(fhat$eval.points[[1]], fhat$eval.points[[2]], 
					fhat$estimate * scale, level = hts[i] * scale, 
					add = TRUE, drawlabels = drawlabels, labcex = labcex, 
					col = col[i], lwd = lwd, ...)
		}
		if (drawpoints) 
			points(fhat$x[, 1], fhat$x[, 2], col = ptcol, cex = cex, 
				pch = pch)
	}
	else if (disp1 == "i") {
		image(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate, 
			xlab = xlab, ylab = ylab, add = add, ...)
		box()
	}
	else if (disp1 == "f") {
		if (display == "filled.contour2") {
			if (missing(abs.cont)) {
				if (!is.null(fhat$cont)) {
					cont.ind <- rep(FALSE, length(fhat$cont))
					for (j in 1:length(cont)) cont.ind[which(cont[j] == 
									as.numeric(unlist(strsplit(names(fhat$cont), 
												"%"))))] <- TRUE
					if (all(!cont.ind)) 
						hts <- contourLevels(fhat, prob = (100 - 
									cont)/100, approx = approx.cont)
					else hts <- fhat$cont[cont.ind]
				}
				else hts <- contourLevels(fhat, prob = (100 - 
								cont)/100, approx = approx.cont)
			}
			else hts <- abs.cont
			hts <- sort(hts)
			if (missing(col)) 
				col <- c("transparent", rev(heat.colors(length(hts))))
			clev <- c(-0.01 * max(abs(fhat$estimate)), hts, max(c(fhat$estimate, 
						hts)) + 0.01 * max(abs(fhat$estimate)))
			image(fhat$eval.points[[1]], fhat$eval.points[[2]], 
				fhat$estimate, xlab = xlab, ylab = ylab, add = add, 
				col = col[1:(length(hts) + 1)], breaks = clev, 
				...)
			for (i in 1:length(hts)) contour(fhat$eval.points[[1]], 
					fhat$eval.points[[2]], fhat$estimate, level = hts[i], 
					add = TRUE, drawlabels = FALSE, col = col[i + 
							1], lwd = 7)
			if (!missing(lwd)) {
				for (i in 1:length(hts)) {
					if (missing(abs.cont)) 
						scale <- cont[i]/hts[i]
					else scale <- 1
					contour(fhat$eval.points[[1]], fhat$eval.points[[2]], 
						fhat$estimate * scale, level = hts[i] * scale, 
						add = TRUE, drawlabels = drawlabels, col = 1, 
						labcex = labcex, lwd = lwd, ...)
				}
			}
		}
		else filled.contour(fhat$eval.points[[1]], fhat$eval.points[[2]], 
				fhat$estimate, xlab = xlab, ylab = ylab, ...)
	}
	if (disp1 == "p") 
		invisible(plotret)
	else invisible()
}

filled.contour <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
		length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
	ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
	levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
	col = color.palette(length(levels) - 1), plot.title, plot.axes, 
	key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
	axes = TRUE, frame.plot = axes, ...) 
{
	if (missing(z)) {
		if (!missing(x)) {
			if (is.list(x)) {
				z <- x$z
				y <- x$y
				x <- x$x
			}
			else {
				z <- x
				x <- seq.int(0, 1, length.out = nrow(z))
			}
		}
		else stop("no 'z' matrix specified")
	}
	else if (is.list(x)) {
		y <- x$y
		x <- x$x
	}
	if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
		stop("increasing 'x' and 'y' values expected")
	mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
	on.exit(par(par.orig))
	w <- (3 + mar.orig[2L]) * par("csi") * 2.54
	layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
	par(las = las)
	mar <- mar.orig
	mar[4L] <- mar[2L]
	mar[2L] <- 1
	par(mar = mar)
	plot.new()
	plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
		yaxs = "i")
	rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
	if (missing(key.axes)) {
		if (axes) 
			axis(4)
	}
	else key.axes
	box()
	if (!missing(key.title)) 
		key.title
	mar <- mar.orig
	mar[4L] <- 1
	par(mar = mar)
	plot.new()
	plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
	if (!is.matrix(z) || nrow(z) <= 1L || ncol(z) <= 1L) 
		stop("no proper 'z' matrix specified")
	if (!is.double(z)) 
		storage.mode(z) <- "double"
	.Internal(filledcontour(as.double(x), as.double(y), z, as.double(levels), 
			col = col))
	if (missing(plot.axes)) {
		if (axes) {
			title(main = "", xlab = "", ylab = "")
			Axis(x, side = 1)
			Axis(y, side = 2)
		}
	}
	else plot.axes
	if (frame.plot) 
		box()
	if (missing(plot.title)) 
		title(...)
	else plot.title
	invisible()
}

