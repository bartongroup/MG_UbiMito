options(dplyr.summarise.inform = FALSE)

files_R <- list.files("R", pattern="*.R$", full.names=TRUE)
sr_ <- sapply(files_R, source)

files_drake <- list.files("drake", pattern="*.R$", full.names=TRUE)
sd_ <- sapply(files_drake, source)

select <- dplyr::select

for(d in c("data", "shiny_all")) {
  if(!dir.exists(d)) dir.create(d)
}

sesinfo <- drake_plan(
  session_info = sessionInfo()
)

plan <- bind_rows(
  get_biomart,
  get_data,
  process_data,
  get_numbers,
  compare_data,
  make_figures,
  save_tables,
  save_for_shiny,
  sesinfo
)

cfg <- drake_config(plan)