library(shiny)
library(wesanderson)
library(dplyr)
COL <- wes_palette("Cavalcanti", 5)
COL[5] <- '#f03b20'
COL[4] <- '#bd0026'
d <- density(faithful$eruptions, adjust = .25)


branchNbound <- function(m, d, it, epsilon) {
  i <- 1
  repeat {
    names(m)[c(2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)] <-
      names(m)[c(2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)][c(2, 1, 5, 6, 3, 4, 9, 10, 7, 8, 12, 11, 15, 13, 14)]

    m <- m %>% group_by(membership_old) %>% # branch
      mutate(membership_new = ifelse(active_new,
                                     ifelse((membership_old - size) < ind & ind <= (membership_old - size/2),
                                            membership_old - size/2, membership_old), membership_old)) %>%
      mutate(size = ifelse(active_new, size/2, size)) %>%
      mutate(xmin_new = ifelse(membership_new == membership_old & active_new == 1, (xmax_old - xmin_old)/2 + xmin_old, xmin_old)) %>%
      mutate(xmax_new = ifelse(membership_new < membership_old & active_new == 1, (xmax_old - xmin_old)/2 + xmin_old, xmax_old)) %>%
      ungroup() %>%
      group_by(membership_new) %>% # bound
      mutate(L_new = min(d$y[xmin_new < d$x & d$x < xmax_new])) %>%
      mutate(U_new = max(d$y[xmin_new < d$x & d$x < xmax_new])) %>%
      mutate(envelope_new = U_new - L_new) %>% # kill
      mutate(active_newest = ifelse(envelope_new < epsilon, 0, 1)) %>%
      ungroup()

    i <- i + 1
    if (i > it) break
  }
  m
}

function(input, output, session) {
  observe({
    which_button <- c("option1", "option2", "option3")[(input$it %% 3) + 1]
    updateRadioButtons(session, "radio", selected = which_button)

    new_max <- input$maxiter * 3
    updateSliderInput(session, "it", max = new_max)
  })

  m <- reactive({
    maxit <- input$maxiter
    data.frame("ind" = 1:2^maxit,
                    "membership_old" = rep(2^maxit, 2^maxit),
                    "membership_new" = rep(2^maxit, 2^maxit),
                    "size" = rep(2^maxit, 2^maxit),
                    "xmin_old" = rep(min(d$x), 2^maxit),
                    "xmax_old" = rep(max(d$x), 2^maxit),
                    "xmin_new" = rep(min(d$x), 2^maxit),
                    "xmax_new" = rep(max(d$x), 2^maxit),
                    "L_old" = rep(min(d$y), 2^maxit),
                    "U_old" = rep(max(d$y), 2^maxit),
                    "L_new" = rep(min(d$y), 2^maxit),
                    "U_new" = rep(max(d$y), 2^maxit),
                    "envelope_old" = rep(max(d$y) - min(d$y), 2^maxit),
                    "envelope_new" = rep(max(d$y) - min(d$y), 2^maxit),
                    "active_old" = rep(1, 2^maxit),
                    "active_new" = rep(1, 2^maxit),
                    "active_newest" = rep(1, 2^maxit))
  })


  output$bbplot <- renderPlot({
    plot(range(d$x), range(d$y), type = "n", xlab = "x", ylab = "f(x)",
         bty = "n", xaxp = c(1, 5.5, 10))

    m2 <- branchNbound(m(), d, it = ceiling(input$it/3), epsilon = input$eps)
    segs <- unique(m2$membership_old)

    for (i in segs) { # bound
      seg_row <- m2 %>% filter(membership_old == i) %>% slice(1)
      lines(c(seg_row$xmin_old, seg_row$xmax_old), rep(seg_row$L_old, 2), col = COL[2], lwd = 3)
      lines(c(seg_row$xmin_old, seg_row$xmax_old), rep(seg_row$U_old, 2), col = COL[2], lwd = 3)
      if(seg_row$active_old) {
        lines(rep(seg_row$xmax_old, 2), c(-.2, 1.05 * max(d$y)), lty = 2, col = "lightgrey", lwd = 2) # (old branch)
      }
      if (seg_row$active_old + seg_row$active_new == 0) { # (old kill)
        seg_row_end <- m2 %>% filter(membership_old == i) %>% slice(n())
        polygon(c(seg_row$xmin_old, seg_row_end$xmax_old, seg_row_end$xmax_old, seg_row$xmin_old),
                c(seg_row$L_old, seg_row$L_old, seg_row$U_old, seg_row$U_old),
                col = COL[4], border = NA)
      }
      if (input$it %in% seq(2, 24, 3)) { # new kill
        if (seg_row$active_old + seg_row$active_new == 1) {
          seg_row_end <- m2 %>% filter(membership_old == i) %>% slice(n())
          polygon(c(seg_row$xmin_old, seg_row_end$xmax_old, seg_row_end$xmax_old, seg_row$xmin_old),
                  c(seg_row$L_old, seg_row$L_old, seg_row$U_old, seg_row$U_old),
                  col = COL[5], border = NA)
        }
      }
      if (input$it %in% seq(3, 24, 3)) { # recent kill
        if (seg_row$active_old + seg_row$active_new == 1) {
          seg_row_end <- m2 %>% filter(membership_old == i) %>% slice(n())
          polygon(c(seg_row$xmin_old, seg_row_end$xmax_old, seg_row_end$xmax_old, seg_row$xmin_old),
                  c(seg_row$L_old, seg_row$L_old, seg_row$U_old, seg_row$U_old),
                  col = COL[4], border = NA)
        }
      }
    }
    lines(d, col = COL[1], lwd = 3)

    if (input$it %in% seq(3, 24, 3)) { # branch
      seg_split <- m2$membership_new[!m2$membership_new == m2$membership_old]
      for(i in seg_split) {
        seg_row <- m2 %>% filter(membership_new == i) %>% slice(1)
        lines(rep(seg_row$xmax_new, 2), c(-.2, 1.05 * max(d$y)), lty = 2, col = "darkgrey", lwd = 2)
      }
    }

  })
  }
