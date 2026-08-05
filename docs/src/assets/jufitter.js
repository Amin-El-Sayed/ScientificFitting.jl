(function () {
  const allowedThemes = new Set(["documenter-light", "documenter-dark"]);
  const themeStorageKey = "documenter-theme";
  const plotStyleStorageKey = "jufitter-plot-style";
  const themeLabels = {
    "documenter-light": "light",
    "documenter-dark": "dark",
  };
  const plotStyles = ["lab", "screen", "article"];
  const legacyLabStyles = new Set(["workbench"]);
  const legacyScreenStyles = new Set(["modern", "showcase"]);
  const defaultPlotStyle = "screen";

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
      if (legacyLabStyles.has(stored)) {
        window.localStorage.setItem(plotStyleStorageKey, "lab");
        return "lab";
      }
      if (legacyScreenStyles.has(stored)) {
        window.localStorage.setItem(plotStyleStorageKey, "screen");
        return "screen";
      }
    } catch (_) {
      // Storage is optional; the selector should still work for this page.
    }
    return defaultPlotStyle;
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
      applyPlotStyle(select.value);
    });

    label.appendChild(select);
    anchor.appendChild(label);
  }

  function applyPlotStyle(style) {
    const selected = plotStyles.includes(style) ? style : defaultPlotStyle;
    document.documentElement.dataset.jufitterPlotStyle = selected;
    document.documentElement.classList.add("jufitter-js-plot-style-ready");

    document.querySelectorAll(".jufitter-plot-style-select").forEach((select) => {
      select.value = selected;
    });

    const groups = new Map();
    document.querySelectorAll("[data-jufitter-plot-group][data-jufitter-plot-style]").forEach((node) => {
      const group = node.getAttribute("data-jufitter-plot-group");
      if (!groups.has(group)) groups.set(group, []);
      groups.get(group).push(node);
    });

    groups.forEach((nodes) => {
      const selectedNodes = nodes.filter((node) => node.getAttribute("data-jufitter-plot-style") === selected);
      const fallbackNodes = selectedNodes.length
        ? selectedNodes
        : nodes.filter((node) => node.getAttribute("data-jufitter-plot-style") === defaultPlotStyle);
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
    buildThemePicker(anchor);
    buildPlotStylePicker(anchor);
    applyPlotStyle(selectedPlotStyle());
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initialize);
  } else {
    initialize();
  }
})();
