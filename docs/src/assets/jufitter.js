(function () {
  const allowedThemes = new Set(["documenter-light", "documenter-dark"]);
  const themeStorageKey = "documenter-theme";
  const plotStyleStorageKey = "jufitter-plot-style";
  const plotPanelStorageKey = "jufitter-plot-panel";
  const themeLabels = {
    "documenter-light": "light",
    "documenter-dark": "dark",
  };
  const plotStyles = ["sans", "tex"];
  const plotPanels = ["show", "hide"];
  const legacyPlotStyles = new Map([
    ["analysis", "sans"],
    ["presentation", "sans"],
    ["screen", "sans"],
    ["lab", "sans"],
    ["workbench", "sans"],
    ["modern", "sans"],
    ["clean", "sans"],
    ["minimal", "sans"],
    ["showcase", "sans"],
    ["article", "tex"],
  ]);
  const defaultPlotStyle = "sans";
  const defaultPlotPanel = "show";

  try {
    const selectedTheme = window.localStorage && window.localStorage.getItem(themeStorageKey);
    if (selectedTheme && selectedTheme.startsWith("catppuccin")) {
      const fallbackTheme = selectedTheme === "catppuccin-latte" ? "documenter-light" : "documenter-dark";
      window.localStorage.setItem(themeStorageKey, fallbackTheme);
      window.location.reload();
      return;
    }
  } catch (_) {
    // Theme cleanup is cosmetic; documentation must still work if storage is blocked.
  }

  function selectedTheme() {
    try {
      const stored = window.localStorage && window.localStorage.getItem(themeStorageKey);
      if (allowedThemes.has(stored)) return stored;
    } catch (_) {
      // Fall through to the DOM class when storage is blocked.
    }
    if (document.documentElement.classList.contains("theme--documenter-dark")) {
      return "documenter-dark";
    }
    return "documenter-light";
  }

  function cleanThemePicker() {
    const picker = document.getElementById("documenter-themepicker");
    if (!picker) return;

    Array.from(picker.options).forEach((option) => {
      if (!allowedThemes.has(option.value)) option.remove();
    });
  }

  function removeUnusedDocumenterToolbarButtons() {
    const selectors = [
      ".docs-settings-button",
      ".docs-article-toggle-button",
      "#documenter-settings-button",
      "#documenter-article-toggle-button",
      "button[aria-label='Settings']",
      "a[aria-label='Settings']",
    ];

    selectors.forEach((selector) => {
      document.querySelectorAll(selector).forEach((node) => node.remove());
    });
  }

  function buildThemePicker(anchor) {
    if (document.querySelector(".jufitter-theme-control")) return;

    const label = document.createElement("label");
    label.className = "jufitter-toolbar-control jufitter-theme-control";
    label.textContent = "Theme";

    const select = document.createElement("select");
    select.className = "jufitter-toolbar-select jufitter-theme-select";
    for (const theme of allowedThemes) {
      const option = document.createElement("option");
      option.value = theme;
      option.textContent = themeLabels[theme];
      select.appendChild(option);
    }

    select.value = selectedTheme();
    select.addEventListener("change", () => {
      try {
        window.localStorage && window.localStorage.setItem(themeStorageKey, select.value);
      } catch (_) {
        // Reload still lets Documenter apply its default theme.
      }
      window.location.reload();
    });

    label.appendChild(select);
    anchor.appendChild(label);
  }

  function selectedPlotStyle() {
    try {
      const stored = window.localStorage && window.localStorage.getItem(plotStyleStorageKey);
      if (plotStyles.includes(stored)) return stored;
      if (legacyPlotStyles.has(stored)) {
        const migrated = legacyPlotStyles.get(stored);
        window.localStorage.setItem(plotStyleStorageKey, migrated);
        return migrated;
      }
    } catch (_) {
      // Storage is optional; the selector should still work for this page.
    }
    return defaultPlotStyle;
  }

  function selectedPlotPanel() {
    try {
      const stored = window.localStorage && window.localStorage.getItem(plotPanelStorageKey);
      if (plotPanels.includes(stored)) return stored;

      // Preserve the composition selected by the former three-role picker.
      const legacyStyle = window.localStorage && window.localStorage.getItem(plotStyleStorageKey);
      if (legacyStyle === "presentation" || legacyStyle === "article") return "hide";
    } catch (_) {
      // Storage is optional; use the self-contained default figure.
    }
    return defaultPlotPanel;
  }

  function toolbarAnchor() {
    const anchor =
      document.querySelector("#documenter .docs-main .docs-navbar .docs-right") ||
      document.querySelector("#documenter .docs-main .docs-navbar") ||
      document.body;
    return anchor;
  }

  function buildPlotStylePicker(anchor) {
    if (document.querySelector(".jufitter-plot-style-control")) return;

    const label = document.createElement("label");
    label.className = "jufitter-toolbar-control jufitter-plot-style-control";
    label.textContent = "Plot style";

    const select = document.createElement("select");
    select.className = "jufitter-toolbar-select jufitter-plot-style-select";
    for (const style of plotStyles) {
      const option = document.createElement("option");
      option.value = style;
      option.textContent = style;
      select.appendChild(option);
    }

    select.value = selectedPlotStyle();
    select.addEventListener("change", () => {
      try {
        window.localStorage && window.localStorage.setItem(plotStyleStorageKey, select.value);
      } catch (_) {
        // Keep the page interactive even when storage is blocked.
      }
      applyPlotSelection(select.value, selectedPlotPanel());
    });

    label.appendChild(select);
    anchor.appendChild(label);
  }

  function buildPlotPanelPicker(anchor) {
    if (document.querySelector(".jufitter-plot-panel-control")) return;

    const label = document.createElement("label");
    label.className = "jufitter-toolbar-control jufitter-plot-panel-control";
    label.textContent = "Result panel";

    const select = document.createElement("select");
    select.className = "jufitter-toolbar-select jufitter-plot-panel-select";
    for (const panel of plotPanels) {
      const option = document.createElement("option");
      option.value = panel;
      option.textContent = panel;
      select.appendChild(option);
    }

    select.value = selectedPlotPanel();
    select.addEventListener("change", () => {
      try {
        window.localStorage && window.localStorage.setItem(plotPanelStorageKey, select.value);
      } catch (_) {
        // Keep the page interactive even when storage is blocked.
      }
      applyPlotSelection(selectedPlotStyle(), select.value);
    });

    label.appendChild(select);
    anchor.appendChild(label);
  }

  function applyPlotSelection(style, panel) {
    const selected = plotStyles.includes(style) ? style : defaultPlotStyle;
    const selectedPanel = plotPanels.includes(panel) ? panel : defaultPlotPanel;
    document.documentElement.dataset.jufitterPlotStyle = selected;
    document.documentElement.dataset.jufitterPlotPanel = selectedPanel;
    document.documentElement.classList.add("jufitter-js-plot-style-ready");

    document.querySelectorAll(".jufitter-plot-style-select").forEach((select) => {
      select.value = selected;
    });
    document.querySelectorAll(".jufitter-plot-panel-select").forEach((select) => {
      select.value = selectedPanel;
    });

    const groups = new Map();
    document.querySelectorAll("[data-jufitter-plot-group][data-jufitter-plot-style]").forEach((node) => {
      const group = node.getAttribute("data-jufitter-plot-group");
      if (!groups.has(group)) groups.set(group, []);
      groups.get(group).push(node);
    });

    groups.forEach((nodes) => {
      const selectedNodes = nodes.filter((node) => {
        const nodePanel = node.getAttribute("data-jufitter-plot-panel");
        return node.getAttribute("data-jufitter-plot-style") === selected &&
          (!nodePanel || nodePanel === selectedPanel);
      });
      const fallbackNodes = selectedNodes.length
        ? selectedNodes
        : nodes.filter((node) => {
            const nodePanel = node.getAttribute("data-jufitter-plot-panel");
            return node.getAttribute("data-jufitter-plot-style") === defaultPlotStyle &&
              (!nodePanel || nodePanel === defaultPlotPanel);
          });
      const visibleNodes = fallbackNodes.length ? fallbackNodes : [nodes[0]];
      const visibleSet = new Set(visibleNodes);

      nodes.forEach((node) => {
        node.classList.toggle("jufitter-style-visible", visibleSet.has(node));
      });
    });
  }

  function initialize() {
    cleanThemePicker();
    removeUnusedDocumenterToolbarButtons();
    const anchor = toolbarAnchor();
    const panel = selectedPlotPanel();
    const style = selectedPlotStyle();
    buildThemePicker(anchor);
    buildPlotStylePicker(anchor);
    buildPlotPanelPicker(anchor);
    applyPlotSelection(style, panel);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initialize);
  } else {
    initialize();
  }
})();
