(function () {
  const allowedThemes = new Set(["auto", "documenter-light", "documenter-dark"]);
  const storageKey = "documenter-theme";

  try {
    const selectedTheme = window.localStorage && window.localStorage.getItem(storageKey);
    if (selectedTheme && selectedTheme.startsWith("catppuccin")) {
      const fallbackTheme = selectedTheme === "catppuccin-latte" ? "documenter-light" : "documenter-dark";
      window.localStorage.setItem(storageKey, fallbackTheme);
      window.location.reload();
      return;
    }
  } catch (_) {
    // Theme cleanup is cosmetic; documentation must still work if storage is blocked.
  }

  function cleanThemePicker() {
    const picker = document.getElementById("documenter-themepicker");
    if (!picker) return;

    Array.from(picker.options).forEach((option) => {
      if (!allowedThemes.has(option.value)) option.remove();
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", cleanThemePicker);
  } else {
    cleanThemePicker();
  }
})();
