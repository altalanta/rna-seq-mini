# 🚀 Quick Start Reference Card

## First-Time Setup (5 minutes)

```bash
# 1. Clone repository
git clone https://github.com/altalanta/nhanes-bmi-bodyfat.git
cd nhanes-bmi-bodyfat

# 2. Check system health
make health-check

# 3. Launch tutorial (if needed)
make tutorial

# 4. Configure analysis
make config-wizard

# 5. Run analysis
make parallel-pipeline
```

## 📊 Essential Commands

### Analysis Pipeline
```bash
make parallel-pipeline    # Complete analysis (recommended)
make all                 # Traditional sequential pipeline
make fetch              # Download data only
make cleandata          # Process data only
make analysis           # Run analysis only
make viz               # Generate visualizations
make report            # Create HTML report
```

### Data Management
```bash
make data-health       # Check data integrity
make data-updates      # Check for new releases
make data-manifest     # Generate reproducibility manifest
```

### Interactive Tools
```bash
make tutorial             # Interactive learning tutorial
make tutorial-troubleshooting  # Troubleshooting guide
make config-wizard        # Visual configuration interface
make demo                # Performance demonstration
```

### Development & Testing
```bash
make test               # Run test suite
make quality           # Code quality checks
make clean             # Remove output files
make cleanall          # Remove everything
```

## 📁 Key Files & Directories

### Configuration
- `config/config.yml` - Analysis parameters
- `README.md` - Comprehensive documentation

### Data & Results
- `data/raw/` - NHANES source files
- `data/derived/` - Processed datasets
- `outputs/tables/` - Statistical results (CSV)
- `outputs/figures/` - Visualizations (PNG/PDF)
- `outputs/report/` - HTML documentation

### System Files
- `Makefile` - Build automation
- `R/` - R package source code
- `tests/` - Test suite
- `docs/` - Documentation

## 🔧 Configuration Options

### Age Range
```yaml
analysis:
  age_range: [20, 59]  # Default: working adults
```

### Survey Design
```yaml
analysis:
  survey_weights_col: "WTMEC2YR"  # MEC examination weights
  strata_col: "SDMVSTRA"         # Survey strata
  psu_col: "SDMVPSU"            # Primary sampling units
```

### Output Directories
```yaml
outputs:
  tables_dir: "outputs/tables"
  figures_dir: "outputs/figures"
  logs_dir: "outputs/logs"
  report_dir: "outputs/report"
```

## 📈 Expected Results

### Correlation Results
- **Overall BMI-body fat correlation**: ~0.914
- **Male correlation**: ~0.917
- **Female correlation**: ~0.954
- **Sample size**: ~2,240 adults

### Output Files
- `corr_bmi_bodyfat_overall_and_by_sex.csv` - Correlation coefficients
- `bodyfat_by_bmi_class_by_sex.csv` - BMI class statistics
- `bmi_vs_bodyfat_plot_sex_facets.png` - Main visualization
- `report.html` - Complete analysis report

## 🚨 Quick Troubleshooting

### "Package not found"
```bash
R -e "install.packages(c('dplyr', 'ggplot2', 'survey', 'foreign'))"
```

### "Directory not found"
```bash
mkdir -p data/raw data/derived outputs/{tables,figures,logs,report}
```

### "Configuration error"
```bash
make config-wizard  # Creates valid config file
```

## 🎓 Learning Resources

### Interactive Tutorial
```bash
make tutorial  # Complete step-by-step guide
```

### Troubleshooting Guide
```bash
make tutorial-troubleshooting  # Issue resolution
```

### Performance Demo
```bash
make demo  # See parallel processing in action
```

## 📞 Getting Help

- **GitHub Issues**: https://github.com/altalanta/nhanes-bmi-bodyfat/issues
- **Email Support**: analysis@nhanes-bmi.org
- **Documentation**: outputs/report/report.html

---

**Print this card and keep it handy!** 📋✨


