my_theme =
	theme_bw() +
	theme(
		panel.border = element_blank(),
		axis.line = element_line(),
		panel.grid.major = element_line(size = 0.2),
		panel.grid.minor = element_line(size = 0.1),
		text = element_text(size=12),
		legend.position="bottom",
		aspect.ratio=1,
		strip.background = element_blank(),
		axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
		axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
		axis.text.x = element_text(angle = 90, hjust = 1)
	)

stats =
	readRDS("dev/stats_pasilla_ttBulk.rds") %>%
	mutate(method = "ttBulk") %>%
	bind_rows(
		readRDS("dev/stats_pasilla_standard.rds") %>%
		mutate(method = "standard")
	) %>%
	mutate(elapsed = map(time, ~ as.numeric(.x$toc - .x$tic))) %>%
	unnest(elapsed) %>%
	select(-time)

stats %>%
	pivot_longer(names_to = ".variable", values_to = ".value", cols = c(lines, assignments, elapsed)) %>%
	mutate(.variable = .variable %>% factor(levels = c("lines", "assignments", "elapsed"))) %>%
	ggplot(aes(x = step, y = .value, color = method, group = method)) +
	geom_line() +
	geom_point() +
	facet_wrap(~.variable, scales = "free") +
	my_theme
