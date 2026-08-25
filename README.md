# KASPtools

`KASPtools` — R-пакет для обработки KASP/KASP-like генотипирования с двух типов выгрузок:

- **DT96 RDML/XML** — импорт всех каналов и всех точек регистрации, включая конечные измерения при 37 °C;
- **QuantStudio XLSX** — импорт `Results` + `Multicomponent Data` для выбранного SNP assay.

Оба импорта формируют единый объект `kasp_run`, после чего используются общие функции визуализации, plateau correction, кластеризации и экспорта.

## Основные функции

- `import_kasp_xml()` — импорт DT96 RDML/XML;
- `import_kasp_quantstudio()` — импорт QuantStudio XLSX;
- `plot_kasp_kinetics()` — графики флуоресценции по всем точкам регистрации;
- `kasp_cartesian()` / `analyse_kasp()` — allelic discrimination, эмпирические кластеры и confidence;
- `kasp_validation_defaults()` — текущие параметры алгоритма валидации;
- `export_kasp_results()` — Excel + опциональный интерактивный HTML.

## Установка с GitHub

После размещения этого каталога в GitHub-репозитории:

```r
install.packages("remotes")
remotes::install_github("USERNAME/KASPtools")
```

или:

```r
install.packages("pak")
pak::pak("USERNAME/KASPtools")
```

Перед публикацией замените `USERNAME` на имя GitHub-пользователя/организации и исправьте `Authors@R` в `DESCRIPTION`.

## DT96: базовый workflow

```r
library(KASPtools)

kasp <- import_kasp_xml(
  file = "RDML_export.xml",
  protocol_to_keep = 1,
  protocol_name = "Rw24"
)

plot_kasp_kinetics(kasp)
```

После выбора plateau window:

```r
result <- kasp_cartesian(
  kasp = kasp,
  plateau_start = 15,
  plateau_end = 25,
  allele_1_control_wells = "D1",
  allele_2_control_wells = "A1",
  heterozygote_control_wells = "C1",
  ntc_wells = "H3",
  allele_1_name = "S",
  allele_2_name = "RW",
  normalization = "delta"
)

result$plot
```

## Режимы нормализации

```r
normalization = "delta"
```

использует:

`endpoint - plateau`

```r
normalization = "relative_delta"
```

использует:

`(endpoint - plateau) / plateau`

```r
normalization = "rox"
```

использует:

`(endpoint - plateau) / ROX_endpoint`

В режиме `rox` канал ROX должен присутствовать в исходной выгрузке.

## QuantStudio

```r
kasp_qs <- import_kasp_quantstudio(
  file = "QuantStudio_export.xlsx",
  assay_to_keep = "SNP Assay 1",
  protocol_name = "Rf1"
)

plot_kasp_kinetics(kasp_qs)

result_qs <- kasp_cartesian(
  kasp = kasp_qs,
  plateau_start = 10,
  plateau_end = 22,
  allele_1_name = "S",
  allele_2_name = "R",
  allele_2_dye = "VIC",
  normalization = "rox"
)
```

Если QuantStudio экспорт содержит `Task = PC_ALLELE_1`, `PC_ALLELE_2`, `PC_ALLELE_BOTH` и `NTC`, контрольные лунки повторно указывать не нужно.

## Экспорт

```r
export_kasp_results(
  result,
  xlsx_output = "KASP_calls.xlsx",
  html_output = "KASP_plot.html",
  include_coordinates = TRUE
)
```

## Настройка confidence

```r
kasp_validation_defaults()
```

Отдельные параметры можно изменить при вызове:

```r
result <- kasp_cartesian(
  kasp = kasp,
  plateau_start = 15,
  plateau_end = 25,
  allele_1_control_wells = "D1",
  allele_2_control_wells = "A1",
  heterozygote_control_wells = "C1",
  ntc_wells = "H3",
  normalization = "delta",
  validation_settings = list(
    min_cluster_separation = 0.12,
    local_neighbours = 4
  )
)
```

## Размещение в GitHub

В корне пакета уже находится GitHub Actions workflow `.github/workflows/R-CMD-check.yaml`. После push GitHub автоматически установит зависимости и выполнит `R CMD check` на Linux.

Минимальная последовательность команд:

```bash
git init
git add .
git commit -m "Initial KASPtools package"
git branch -M main
git remote add origin git@github.com:USERNAME/KASPtools.git
git push -u origin main
```

## License

MIT. Перед публичным размещением при необходимости измените правообладателя в `LICENSE` и `Authors@R` в `DESCRIPTION`.
