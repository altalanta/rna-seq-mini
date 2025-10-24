# Troubleshooting Decision Tree Flowchart

A visual decision tree to help you quickly diagnose and resolve common issues with the NHANES BMI Body Fat Analysis platform.

```
🔍 START: Experiencing an issue?

├── 📦 PACKAGE/INSTALLATION ISSUES
│   ├── "Package not found" errors?
│   │   ├── ✅ QUICK FIX: Install missing packages
│   │   │   Command: R -e "install.packages(c('dplyr', 'ggplot2', 'survey', 'foreign'))"
│   │   │   └── Continue to "Test installation"
│   │   │
│   │   └── "Installation failed" errors?
│   │       ├── Check R version compatibility
│   │       │   Command: R -e "sessionInfo()"
│   │       │   └── Update R if needed
│   │       │
│   │       ├── Check system dependencies
│   │       │   Linux: sudo apt install libxml2-dev libssl-dev
│   │       │   macOS: brew install libxml2 openssl
│   │       │   └── Windows: Install Rtools
│   │       │
│   │       └── Try alternative installation method
│   │           └── Use: make config-wizard
│   │
│   └── Test installation
│       ├── Run: make health-check
│       ├── Expected: "All systems operational!"
│       └── ✅ If passed → Issue resolved
│
├── 📁 FILE/DIRECTORY ISSUES
│   ├── "Directory not found" errors?
│   │   ├── ✅ QUICK FIX: Create required directories
│   │   │   Command: mkdir -p data/raw data/derived outputs/{tables,figures,logs,report}
│   │   │   └── Continue to "Test directories"
│   │   │
│   │   └── "Permission denied" errors?
│   │       ├── Check file permissions
│   │       │   Command: ls -la
│   │       │   └── Fix: chmod 755 data outputs config
│   │       │
│   │       └── Run in user directory
│   │           └── mkdir -p ~/nhanes-analysis && cd ~/nhanes-analysis
│   │
│   └── Test directories
│       ├── Run: make health-check
│       └── ✅ If passed → Issue resolved
│
├── ⚙️ CONFIGURATION ISSUES
│   ├── "Configuration file not found"?
│   │   ├── ✅ QUICK FIX: Use configuration wizard
│   │   │   Command: make config-wizard
│   │   │   └── Creates config/config.yml automatically
│   │   │
│   │   └── "Invalid YAML" errors?
│   │       ├── Check YAML syntax
│   │       │   Use spaces, not tabs
│   │       │   Proper indentation required
│   │       │   └── Validate: R -e "yaml::read_yaml('config/config.yml')"
│   │       │
│   │       └── Manual creation
│   │           └── Copy from config/config.yml.example
│   │
│   └── Test configuration
│       ├── Run: make health-check
│       └── ✅ If passed → Issue resolved
│
├── 🌐 NETWORK/DOWNLOAD ISSUES
│   ├── "Data download failed"?
│   │   ├── Check internet connection
│   │   │   Command: ping -c 3 google.com
│   │   │   └── If failed → Fix network connection
│   │   │
│   │   ├── Check available disk space
│   │   │   Command: df -h
│   │   │   └── Need ~100MB free space
│   │   │
│   │   ├── Retry download
│   │   │   Command: make fetch
│   │   │   └── Check logs for specific errors
│   │   │
│   │   └── Manual download alternative
│   │       └── Download files manually from CDC website
│   │
│   └── Test download
│       ├── Run: make fetch
│       └── ✅ If completed → Issue resolved
│
├── 📊 ANALYSIS ISSUES
│   ├── "Analysis failed to run"?
│   │   ├── Check system requirements
│   │   │   Memory: 4GB+ recommended
│   │   │   CPU: Multi-core for parallel processing
│   │   │   └── Disk space: 1GB+ available
│   │   │
│   │   ├── Check data availability
│   │   │   Command: ls data/raw/
│   │   │   └── Should contain .XPT files
│   │   │
│   │   ├── Try parallel pipeline
│   │   │   Command: make parallel-pipeline
│   │   │   └── More robust than sequential
│   │   │
│   │   └── Check configuration
│   │       └── Verify age range and parameters
│   │
│   └── Test analysis
│       ├── Run: make parallel-pipeline
│       └── ✅ If completed → Issue resolved
│
├── 💾 DATA INTEGRITY ISSUES
│   ├── "Hash mismatch" errors?
│   │   ├── ✅ QUICK FIX: Re-download data
│   │   │   Command: make fetch
│   │   │   └── Verifies file integrity automatically
│   │   │
│   │   └── "File corrupted" errors?
│   │       ├── Check file size
│   │       │   Command: ls -lh data/raw/
│   │       │   └── Compare with expected sizes
│   │       │
│   │       └── Manual verification
│   │           └── Re-download specific files
│   │
│   └── Test integrity
│       ├── Run: make data-integrity
│       └── ✅ If passed → Issue resolved
│
├── ⚡ PERFORMANCE ISSUES
│   ├── "Analysis too slow"?
│   │   ├── ✅ QUICK FIX: Use parallel processing
│   │   │   Command: make parallel-pipeline
│   │   │   └── 3-5x faster than sequential
│   │   │
│   │   ├── Check system resources
│   │   │   Command: top -o cpu (Linux/macOS)
│   │   │   └── Task Manager (Windows)
│   │   │
│   │   ├── Memory issues?
│   │   │   Command: R -e "memory.limit()"
│   │   │   └── Increase if needed
│   │   │
│   │   └── Hardware recommendations
│   │       └── SSD storage, 8GB+ RAM, multi-core CPU
│   │
│   └── Test performance
│       ├── Run: make demo
│       └── ✅ If improved → Issue resolved
│
├── 🎓 TUTORIAL ISSUES
│   ├── "Tutorial won't launch"?
│   │   ├── ✅ QUICK FIX: Check browser compatibility
│   │   │   Works best in Chrome, Firefox, Safari
│   │   │   Try different browser if needed
│   │   │
│   │   ├── Check learnr package
│   │   │   Command: R -e "library(learnr)"
│   │   │   └── Install if missing
│   │   │
│   │   └── Port conflicts?
│   │       └── Change port: R -e "learnr::run_tutorial('tutorials/getting_started.Rmd', port=3839)"
│   │
│   └── Test tutorial
│       ├── Run: make tutorial
│       └── ✅ If launches → Issue resolved
│
├── ⚙️ CONFIGURATION WIZARD ISSUES
│   ├── "Wizard won't launch"?
│   │   ├── ✅ QUICK FIX: Check Shiny package
│   │   │   Command: R -e "library(shiny)"
│   │   │   └── Install if missing
│   │   │
│   │   ├── Port conflicts?
│   │   │   Command: lsof -i :3838 (Linux/macOS)
│   │   │   └── netstat -ano | findstr :3838 (Windows)
│   │   │
│   │   └── Browser issues?
│   │       └── Try: R -e "shiny::runApp('app.R', launch.browser=FALSE)"
│   │
│   └── Test wizard
│       ├── Run: make config-wizard
│       └── ✅ If launches → Issue resolved
│
└── 🆘 ESCALATION: Need More Help?

    ├── Self-service resources
    │   ├── Interactive tutorial: make tutorial
    │   ├── Troubleshooting guide: make tutorial-troubleshooting
    │   ├── Health check: make health-check
    │   └── Error logs: outputs/logs/
    │
    ├── Community support
    │   ├── GitHub Issues: https://github.com/altalanta/nhanes-bmi-bodyfat/issues
    │   ├── GitHub Discussions: https://github.com/altalanta/nhanes-bmi-bodyfat/discussions
    │   └── Documentation: outputs/report/report.html
    │
    └── Professional support
        ├── Training workshops: Contact analysis@nhanes-bmi.org
        ├── Custom development: Specialized implementations
        ├── Code review: Professional assessment
        └── Consulting: Methodological guidance

🎉 ISSUE RESOLVED? Great! Continue with your research.

❌ STILL HAVING ISSUES? Contact support with:
- Error messages and logs
- System information (R version, OS, etc.)
- Steps to reproduce the issue
- Expected vs actual behavior
```

## 📋 Quick Reference Commands

### Diagnostic Commands
```bash
make health-check          # Overall system check
make data-health          # Data integrity check
make data-integrity       # File integrity validation
make data-updates         # Check for data updates
R -e "sessionInfo()"      # R environment info
R -e "memory.size()"      # Memory usage info
```

### Fix Commands
```bash
# Package issues
R -e "install.packages(c('dplyr', 'ggplot2', 'survey', 'foreign'))"

# Directory issues
mkdir -p data/raw data/derived outputs/{tables,figures,logs,report}

# Configuration issues
make config-wizard

# Data issues
make fetch

# Performance issues
make parallel-pipeline

# Tutorial issues
make tutorial

# API issues
make api-test
```

### Monitoring Commands
```bash
# Real-time monitoring
make monitor

# Performance tracking
make performance-tools

# System health
make health-check

# Data quality
make data-health
```

## 🚨 Emergency Procedures

### Critical System Failure
1. **Run health check:** `make health-check`
2. **Check logs:** `cat outputs/logs/analysis_log.txt | tail -20`
3. **Try recovery:** `make clean-cache && make data-registry-init`
4. **Contact support** if issues persist

### Data Loss
1. **Check backups:** `ls backups/`
2. **Restore from backup:** `Rscript deployment/backup_recovery.R restore [backup_name]`
3. **Verify restoration:** `make data-integrity`
4. **Re-run analysis** if needed

### Performance Crisis
1. **Monitor resources:** `top -o cpu` (Linux/macOS)
2. **Reduce workers:** Edit `config/config.yml` to lower worker count
3. **Use sequential:** `make all` instead of `make parallel-pipeline`
4. **Check hardware:** Ensure adequate RAM and CPU

---

**This decision tree provides a systematic approach to resolving any issues you encounter with the NHANES BMI Body Fat Analysis platform.**


