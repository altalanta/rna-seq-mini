# Quick Start Walkthrough Guide

A step-by-step interactive guide to get you up and running with the NHANES BMI Body Fat Analysis platform in under 10 minutes.

## 🚀 Step 1: Environment Setup (2 minutes)

### What you'll do:
- Verify your system meets requirements
- Install necessary R packages
- Create required directory structure

### Interactive Steps:

1. **Open your terminal** and navigate to your desired workspace
2. **Clone the repository:**
   ```bash
   git clone https://github.com/altalanta/nhanes-bmi-bodyfat.git
   cd nhanes-bmi-bodyfat
   ```

3. **Check system health:**
   ```bash
   make health-check
   ```
   ✅ **Expected Result:** "All systems operational!"

4. **If you see any missing packages:**
   ```bash
   # Install missing packages
   R -e "install.packages(c('dplyr', 'ggplot2', 'survey', 'foreign'))"
   ```

## ⚙️ Step 2: Configuration Setup (2 minutes)

### What you'll do:
- Customize analysis parameters for your research
- Set up data directories
- Configure output preferences

### Interactive Steps:

1. **Launch the configuration wizard:**
   ```bash
   make config-wizard
   ```
   *This opens a web interface in your browser*

2. **In the web interface, set your preferences:**
   - **Age Range:** Adjust sliders (default: 20-59 years)
   - **Survey Weights:** Select appropriate column (default: WTMEC2YR)
   - **Output Directories:** Verify paths are correct

3. **Save your configuration:**
   - Click "Save Configuration"
   - The system creates `config/config.yml`

## 📊 Step 3: Run Your First Analysis (3 minutes)

### What you'll do:
- Execute the complete analysis pipeline
- Monitor progress in real-time
- Generate results and visualizations

### Interactive Steps:

1. **Run the parallel processing pipeline:**
   ```bash
   make parallel-pipeline
   ```

2. **Watch the progress:**
   ```
   ✅ Step 0: Data version management and integrity checks...
   ✅ Step 1: Loading NHANES datasets...
   ✅ Step 5: Computing correlations in parallel...
   ✅ Step 9: Exporting results...
   ✅ Pipeline completed in 22.8 seconds
   ```

3. **Verify completion:**
   ```bash
   ls outputs/
   ```
   **Expected files:**
   - `tables/corr_bmi_bodyfat_overall_and_by_sex.csv`
   - `figures/bmi_vs_bodyfat_plot_sex_facets.png`
   - `report/report.html`

## 📈 Step 4: Explore Your Results (2 minutes)

### What you'll do:
- Examine statistical results
- View visualizations
- Review the complete analysis report

### Interactive Steps:

1. **View statistical results:**
   ```bash
   # Open correlation results
   open outputs/tables/corr_bmi_bodyfat_overall_and_by_sex.csv
   ```
   *Look for correlation coefficients and confidence intervals*

2. **View visualizations:**
   ```bash
   # Open the main plot
   open outputs/figures/bmi_vs_bodyfat_plot_sex_facets.png
   ```
   *Observe BMI-body fat relationships by sex*

3. **Read the complete report:**
   ```bash
   # Open HTML report
   open outputs/report/report.html
   ```
   *Review methodology, results, and interactive visualizations*

## 🎯 Step 5: Customize for Your Research (1 minute)

### What you'll do:
- Modify parameters for your specific research question
- Re-run analysis with custom settings
- Export results for your use case

### Interactive Steps:

1. **Modify configuration:**
   ```bash
   make config-wizard  # Re-open to adjust settings
   ```
   *Try changing age range or other parameters*

2. **Re-run analysis:**
   ```bash
   make parallel-pipeline  # Uses your new configuration
   ```

3. **Export for your needs:**
   ```bash
   # Copy results to your research directory
   cp outputs/tables/*.csv ~/my_research/
   cp outputs/figures/*.png ~/my_research/
   ```

## ✅ Success Checklist

- [ ] **Environment Setup:** Health check passed
- [ ] **Configuration:** Custom parameters saved
- [ ] **Analysis Complete:** Pipeline finished successfully
- [ ] **Results Available:** Files created in outputs/
- [ ] **Understanding:** Can interpret correlation results
- [ ] **Customization:** Modified settings for your research

## 🚨 Common Issues & Quick Fixes

### Issue: "Package not found"
**Quick Fix:**
```bash
R -e "install.packages(c('dplyr', 'ggplot2', 'survey', 'foreign'))"
```

### Issue: "Directory not found"
**Quick Fix:**
```bash
mkdir -p data/raw data/derived outputs/{tables,figures,logs,report}
```

### Issue: "Configuration file not found"
**Quick Fix:**
```bash
make config-wizard  # Creates the config file for you
```

## 🎓 Next Steps for Learning

1. **Complete the interactive tutorial:**
   ```bash
   make tutorial
   ```

2. **Explore the troubleshooting guide:**
   ```bash
   make tutorial-troubleshooting
   ```

3. **Try advanced features:**
   ```bash
   make demo  # See parallel processing in action
   make api-launch  # Access results via API
   ```

## 📞 Getting Help

- **Interactive Tutorial:** `make tutorial`
- **Troubleshooting Guide:** `make tutorial-troubleshooting`
- **Configuration Help:** `make config-wizard`
- **System Health:** `make health-check`
- **Community Support:** https://github.com/altalanta/nhanes-bmi-bodyfat/issues

---

**🎉 Congratulations!** You've successfully completed your first NHANES BMI Body Fat analysis. You're now ready to conduct professional epidemiological research with confidence!

**Time to completion:** ~10 minutes
**Skill level achieved:** Basic analysis and customization
**Next milestone:** Advanced statistical modeling and custom research questions


